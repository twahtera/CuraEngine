/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <cmath> // std::ceil

#include "skin.h"
#include "utils/math.h"
#include "utils/polygonUtils.h"

#define MIN_AREA_SIZE (0.4 * 0.4) 

namespace cura 
{

SkinInfillAreaComputation::SkinInfillAreaComputation(int layer_nr, SliceMeshStorage& mesh, int bottom_layer_count, int top_layer_count, int wall_line_count, const int innermost_wall_line_width, int infill_skin_overlap, int wall_line_width_x, int skin_inset_count, bool no_small_gaps_heuristic, bool process_infill)
: layer_nr(layer_nr)
, mesh(mesh)
, bottom_layer_count(bottom_layer_count)
, top_layer_count(top_layer_count)
, wall_line_count(wall_line_count)
, innermost_wall_line_width(innermost_wall_line_width)
, infill_skin_overlap(infill_skin_overlap)
, wall_line_width_x(wall_line_width_x)
, skin_inset_count(skin_inset_count)
, no_small_gaps_heuristic(no_small_gaps_heuristic)
, process_infill(process_infill)
{
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
Polygons SkinInfillAreaComputation::getInsidePolygons(const SliceLayerPart& part_here, const SliceLayer& layer2)
{
    Polygons result;
    for (const SliceLayerPart& part2 : layer2.parts)
    {
        if (part_here.boundaryBox.hit(part2.boundaryBox))
        {
            unsigned int wall_idx = std::max(0, std::min(wall_line_count, (int) part2.insets.size()) - 1);
            result.add(part2.insets[wall_idx]);
        }
    }
    return result;
};

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
Polygons SkinInfillAreaComputation::getOutlines(const SliceLayerPart& part_here, int layer2_nr)
{
    Polygons result;
    if (layer2_nr < static_cast<int>(mesh.layers.size()))
    {
        const SliceLayer& layer2 = mesh.layers[layer2_nr];
        for (const SliceLayerPart& part2 : layer2.parts)
        {
            if (part_here.boundaryBox.hit(part2.boundaryBox))
            {
                result.add(part2.print_outline);
            }
        }
    }
    return result;
};

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * generateSkinAreas reads data from mesh.layers.parts[*].insets and writes to mesh.layers[n].parts[*].skin_parts
 * generateSkinInsets only read/writes the skin_parts from the current layer.
 *
 * generateSkins therefore reads (depends on) data from mesh.layers[*].parts[*].insets and writes mesh.layers[n].parts[*].skin_parts
 */
void SkinInfillAreaComputation::generateSkinsAndInfill()
{
    generateSkinAreas();

    SliceLayer* layer = &mesh.layers[layer_nr];
    for (unsigned int part_nr = 0; part_nr < layer->parts.size(); part_nr++)
    {
        SliceLayerPart* part = &layer->parts[part_nr];
        generateSkinInsets(part);
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * generateSkinAreas reads data from mesh.layers[*].parts[*].insets and writes to mesh.layers[n].parts[*].skin_parts
 */
void SkinInfillAreaComputation::generateSkinAreas()
{
    SliceLayer& layer = mesh.layers[layer_nr];
    
    if (bottom_layer_count == 0 && top_layer_count == 0)
    {
        return;
    }

    for (unsigned int part_nr = 0; part_nr < layer.parts.size(); part_nr++)
    {
        SliceLayerPart& part = layer.parts[part_nr];

        if (static_cast<int>(part.insets.size()) < wall_line_count)
        {
            continue; // the last wall is not present, the part should only get inter perimeter gaps, but no skin or infill.
        }
        generateSkinAndInfillAreas(part);
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * generateSkinAreas reads data from mesh.layers[*].parts[*].insets and writes to mesh.layers[n].parts[*].skin_parts
 */
void SkinInfillAreaComputation::generateSkinAndInfillAreas(SliceLayerPart& part)
{
    int min_infill_area = mesh.getSettingInMillimeters("min_infill_area");

    Polygons upskin = part.insets.back().offset(-innermost_wall_line_width / 2);
    // make a copy of the outline which we later intersect and union with the resized skins to ensure the resized skin isn't too large or removed completely.
    Polygons original_outline = Polygons(upskin);
    Polygons downskin = (bottom_layer_count == 0) ? Polygons() : upskin;
    if (top_layer_count == 0)
    {
        upskin = Polygons();
    }

    calculateBottomSkin(part, min_infill_area, downskin);

    calculateTopSkin(part, min_infill_area, upskin);

    applySkinExpansion(original_outline, upskin, downskin);

    // now combine the resized upskin and downskin
    Polygons skin = upskin.unionPolygons(downskin);

    skin.removeSmallAreas(MIN_AREA_SIZE);

    if (process_infill)
    { // process infill when infill density > 0
        // or when other infill meshes want to modify this infill
        generateInfill(part, skin);
    }

    calculateTopMostSkin(part, skin);

    for (PolygonsPart& skin_area_part : skin.splitIntoParts())
    {
        part.skin_parts.emplace_back();
        part.skin_parts.back().outline = skin_area_part;
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
void SkinInfillAreaComputation::calculateBottomSkin(const SliceLayerPart& part, int min_infill_area, Polygons& downskin)
{
    if (static_cast<int>(layer_nr - bottom_layer_count) >= 0 && bottom_layer_count > 0)
    {
        Polygons not_air = getInsidePolygons(part, mesh.layers[layer_nr - bottom_layer_count]);
        if (!no_small_gaps_heuristic)
        {
            for (int downskin_layer_nr = layer_nr - bottom_layer_count + 1; downskin_layer_nr < layer_nr; downskin_layer_nr++)
            {
                not_air = not_air.intersection(getInsidePolygons(part, mesh.layers[downskin_layer_nr]));
            }
        }
        if (min_infill_area > 0)
        {
            not_air.removeSmallAreas(min_infill_area);
        }
        downskin = downskin.difference(not_air); // skin overlaps with the walls
    }
}

void SkinInfillAreaComputation::calculateTopSkin(const SliceLayerPart& part, int min_infill_area, Polygons& upskin)
{
    if (static_cast<int>(layer_nr + top_layer_count) < static_cast<int>(mesh.layers.size()) && top_layer_count > 0)
    {
        Polygons not_air = getInsidePolygons(part, mesh.layers[layer_nr + top_layer_count]);
        if (!no_small_gaps_heuristic)
        {
            for (int upskin_layer_nr = layer_nr + 1; upskin_layer_nr < layer_nr + top_layer_count; upskin_layer_nr++)
            {
                not_air = not_air.intersection(getInsidePolygons(part, mesh.layers[upskin_layer_nr]));
            }
        }
        if (min_infill_area > 0)
        {
            not_air.removeSmallAreas(min_infill_area);
        }
        upskin = upskin.difference(not_air); // skin overlaps with the walls
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
void SkinInfillAreaComputation::applySkinExpansion(const Polygons& original_outline, Polygons& upskin, Polygons& downskin)
{
    coord_t expand_skins_expand_distance = mesh.getSettingInMicrons("expand_skins_expand_distance");
    if (expand_skins_expand_distance <= 0)
    {
        return;
    }

    coord_t pre_shrink = mesh.getSettingInMicrons("min_skin_width_for_expansion") / 2;

    // skin areas are to be enlarged by expand_skins_expand_distance but before they are expanded
    // the skin areas are shrunk by pre_shrink so that very narrow regions of skin
    // (often caused by the model's surface having a steep incline) are removed first

    expand_skins_expand_distance += pre_shrink; // increase the expansion distance to compensate for the shrinkage

    if (mesh.getSettingBoolean("expand_upper_skins"))
    {
        upskin = upskin.offset(-pre_shrink).offset(expand_skins_expand_distance).unionPolygons(upskin).intersection(original_outline);
    }

    if (mesh.getSettingBoolean("expand_lower_skins"))
    {
        downskin = downskin.offset(-pre_shrink).offset(expand_skins_expand_distance).unionPolygons(downskin).intersection(original_outline);
    }
}


/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
void SkinInfillAreaComputation::generateSkinInsets(SliceLayerPart* part)
{
    if (skin_inset_count == 0)
    {
        return;
    }
    
    for (SkinPart& skin_part : part->skin_parts)
    {
        generateSkinInsets(skin_part);
    }
    for (SkinPart& topmost_skin_part : part->topmost_skin_parts)
    {
        generateSkinInsets(topmost_skin_part);
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
void SkinInfillAreaComputation::generateSkinInsets(SkinPart& skin_part)
{
    for (int inset_idx = 0; inset_idx < insetCount; inset_idx++)
    {
        skin_part.insets.push_back(Polygons());
        if (inset_idx == 0)
        {
            skin_part.insets[0] = skin_part.outline.offset(-wall_line_width_x / 2);
        }
        else
        {
            skin_part.insets[inset_idx] = skin_part.insets[inset_idx - 1].offset(-wall_line_width_x);
        }

        // optimize polygons: remove unnecessary verts
        skin_part.insets[inset_idx].simplify();
        if (skin_part.insets[inset_idx].size() < 1)
        {
            skin_part.insets.pop_back();
            break;
        }
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * generateInfill read mesh.layers[n].parts[*].{insets,skin_parts,boundingBox} and write mesh.layers[n].parts[*].infill_area
 */
void SkinInfillAreaComputation::generateInfill(SliceLayerPart& part, const Polygons& skin)
{
    int extra_offset = 0;
    EFillMethod fill_pattern = mesh.getSettingAsFillMethod("infill_pattern");
    if ((fill_pattern == EFillMethod::CONCENTRIC || fill_pattern == EFillMethod::CONCENTRIC_3D)
        && mesh.getSettingBoolean("alternate_extra_perimeter")
        && layer_nr % 2 == 0
        && mesh.getSettingInMicrons("infill_line_distance") > mesh.getSettingInMicrons("infill_line_width") * 2)
    {
        extra_offset = -innermost_wall_line_width;
    }

    if (int(part.insets.size()) < wall_line_count)
    {
        return; // the last wall is not present, the part should only get inter preimeter gaps, but no infill.
    }
    Polygons infill = part.insets.back().offset(extra_offset - innermost_wall_line_width / 2 - infill_skin_overlap);

    infill = infill.difference(skin);
    infill.removeSmallAreas(MIN_AREA_SIZE);

    Polygons final_infill = infill.offset(infill_skin_overlap);

    if (mesh.getSettingBoolean("infill_hollow"))
    {
        part.print_outline = part.print_outline.difference(final_infill);
    }
    else
    {
        part.infill_area = final_infill;
    }
}

/*
 * This function is executed in a parallel region based on layer_nr.
 * When modifying make sure any changes does not introduce data races.
 *
 * this function may only read/write the skin and infill from the *current* layer.
 */
void SkinInfillAreaComputation::calculateTopMostSkin(SliceLayerPart& part, Polygons& skin)
{
    int topmost_skin_layer_count = mesh.getSettingAsCount("topmost_skin_layer_count");

    Polygons topmost_skin;
    if (topmost_skin_layer_count > 0)
    {
        Polygons no_air_above = getOutlines(part, layer_nr + topmost_skin_layer_count);
        if (!no_small_gaps_heuristic)
        {
            for (int layer_nr_above = layer_nr + 1; layer_nr_above < layer_nr + topmost_skin_layer_count; layer_nr_above++)
            {
                Polygons outlines_above = getOutlines(part, layer_nr_above);
                no_air_above.intersection(outlines_above);
            }
        }
        topmost_skin = skin.difference(no_air_above);
        skin = skin.intersection(no_air_above);

        for (PolygonsPart& topmost_skin_area_part : topmost_skin.splitIntoParts())
        {
            part.topmost_skin_parts.emplace_back();
            part.topmost_skin_parts.back().outline = topmost_skin_area_part;
        }
    }
}

void SkinInfillAreaComputation::generateGradualInfill(SliceMeshStorage& mesh, unsigned int gradual_infill_step_height, unsigned int max_infill_steps)
{
    // no early-out for this function; it needs to initialize the [infill_area_per_combine_per_density]
    float layer_skip_count = 8; // skip every so many layers as to ignore small gaps in the model making computation more easy
    if (!mesh.getSettingBoolean("skin_no_small_gaps_heuristic"))
    {
        layer_skip_count = 1;
    }
    unsigned int gradual_infill_step_layer_count = round_divide(gradual_infill_step_height, mesh.getSettingInMicrons("layer_height")); // The difference in layer count between consecutive density infill areas

    // make gradual_infill_step_height divisable by layer_skip_count
    float n_skip_steps_per_gradual_step = std::max(1.0f, std::ceil(gradual_infill_step_layer_count / layer_skip_count)); // only decrease layer_skip_count to make it a divisor of gradual_infill_step_layer_count
    layer_skip_count = gradual_infill_step_layer_count / n_skip_steps_per_gradual_step;


    size_t min_layer = mesh.getSettingAsCount("bottom_layers");
    size_t max_layer = mesh.layers.size() - 1 - mesh.getSettingAsCount("top_layers");

    for (size_t layer_idx = 0; layer_idx < mesh.layers.size(); layer_idx++)
    { // loop also over layers which don't contain infill cause of bottom_ and top_layer to initialize their infill_area_per_combine_per_density
        SliceLayer& layer = mesh.layers[layer_idx];

        for (SliceLayerPart& part : layer.parts)
        {
            assert(part.infill_area_per_combine_per_density.size() == 0 && "infill_area_per_combine_per_density is supposed to be uninitialized");

            const Polygons& infill_area = part.getOwnInfillArea();

            if (infill_area.size() == 0 || layer_idx < min_layer || layer_idx > max_layer)
            { // initialize infill_area_per_combine_per_density empty
                part.infill_area_per_combine_per_density.emplace_back(); // create a new infill_area_per_combine
                part.infill_area_per_combine_per_density.back().emplace_back(); // put empty infill area in the newly constructed infill_area_per_combine
                // note: no need to copy part.infill_area, cause it's the empty vector anyway
                continue;
            }
            Polygons less_dense_infill = infill_area; // one step less dense with each infill_step
            for (unsigned int infill_step = 0; infill_step < max_infill_steps; infill_step++)
            {
                size_t min_layer = layer_idx + infill_step * gradual_infill_step_layer_count + layer_skip_count;
                size_t max_layer = layer_idx + (infill_step + 1) * gradual_infill_step_layer_count;

                for (float upper_layer_idx = min_layer; static_cast<unsigned int>(upper_layer_idx) <= max_layer; upper_layer_idx += layer_skip_count)
                {
                    if (static_cast<unsigned int>(upper_layer_idx) >= mesh.layers.size())
                    {
                        less_dense_infill.clear();
                        break;
                    }
                    SliceLayer& upper_layer = mesh.layers[static_cast<unsigned int>(upper_layer_idx)];
                    Polygons relevent_upper_polygons;
                    for (SliceLayerPart& upper_layer_part : upper_layer.parts)
                    {
                        if (!upper_layer_part.boundaryBox.hit(part.boundaryBox))
                        {
                            continue;
                        }
                        relevent_upper_polygons.add(upper_layer_part.getOwnInfillArea());
                    }
                    less_dense_infill = less_dense_infill.intersection(relevent_upper_polygons);
                }
                if (less_dense_infill.size() == 0)
                {
                    break;
                }
                // add new infill_area_per_combine for the current density
                part.infill_area_per_combine_per_density.emplace_back();
                std::vector<Polygons>& infill_area_per_combine_current_density = part.infill_area_per_combine_per_density.back();
                const Polygons more_dense_infill = infill_area.difference(less_dense_infill);
                infill_area_per_combine_current_density.push_back(more_dense_infill);

                if (less_dense_infill.size() == 0)
                {
                    break;
                }
            }
            part.infill_area_per_combine_per_density.emplace_back();
            std::vector<Polygons>& infill_area_per_combine_current_density = part.infill_area_per_combine_per_density.back();
            infill_area_per_combine_current_density.push_back(infill_area);
            part.infill_area_own = nullptr; // clear infill_area_own, it's not needed any more.
            assert(part.infill_area_per_combine_per_density.size() != 0 && "infill_area_per_combine_per_density is now initialized");
        }
    }
}

void SkinInfillAreaComputation::combineInfillLayers(SliceMeshStorage& mesh, unsigned int amount)
{
    if (mesh.layers.empty() || mesh.layers.size() - 1 < static_cast<size_t>(mesh.getSettingAsCount("top_layers")) || mesh.getSettingAsCount("infill_line_distance") <= 0) //No infill is even generated.
    {
        return;
    }
    if(amount <= 1) //If we must combine 1 layer, nothing needs to be combined. Combining 0 layers is invalid.
    {
        return;
    }
    /* We need to round down the layer index we start at to the nearest
    divisible index. Otherwise we get some parts that have infill at divisible
    layers and some at non-divisible layers. Those layers would then miss each
    other. */
    size_t min_layer = mesh.getSettingAsCount("bottom_layers") + amount - 1;
    min_layer -= min_layer % amount; //Round upwards to the nearest layer divisible by infill_sparse_combine.
    size_t max_layer = mesh.layers.size() - 1 - mesh.getSettingAsCount("top_layers");
    max_layer -= max_layer % amount; //Round downwards to the nearest layer divisible by infill_sparse_combine.
    for(size_t layer_idx = min_layer;layer_idx <= max_layer;layer_idx += amount) //Skip every few layers, but extrude more.
    {
        SliceLayer* layer = &mesh.layers[layer_idx];
        for(unsigned int combine_count_here = 1; combine_count_here < amount; combine_count_here++)
        {
            if(layer_idx < combine_count_here)
            {
                break;
            }

            size_t lower_layer_idx = layer_idx - combine_count_here;
            if (lower_layer_idx < min_layer)
            {
                break;
            }
            SliceLayer* lower_layer = &mesh.layers[lower_layer_idx];
            for (SliceLayerPart& part : layer->parts)
            {
                for (unsigned int density_idx = 0; density_idx < part.infill_area_per_combine_per_density.size(); density_idx++)
                { // go over each density of gradual infill (these density areas overlap!)
                    std::vector<Polygons>& infill_area_per_combine = part.infill_area_per_combine_per_density[density_idx];
                    Polygons result;
                    for (SliceLayerPart& lower_layer_part : lower_layer->parts)
                    {
                        if (part.boundaryBox.hit(lower_layer_part.boundaryBox))
                        {

                            Polygons intersection = infill_area_per_combine[combine_count_here - 1].intersection(lower_layer_part.infill_area).offset(-200).offset(200);
                            result.add(intersection); // add area to be thickened
                            infill_area_per_combine[combine_count_here - 1] = infill_area_per_combine[combine_count_here - 1].difference(intersection); // remove thickened area from less thick layer here
                            if (density_idx < lower_layer_part.infill_area_per_combine_per_density.size())
                            { // only remove from *same density* areas on layer below
                                // If there are no same density areas, then it's ok to print them anyway
                                // Don't remove other density areas
                                unsigned int lower_density_idx = density_idx;
                                std::vector<Polygons>& lower_infill_area_per_combine = lower_layer_part.infill_area_per_combine_per_density[lower_density_idx];
                                lower_infill_area_per_combine[0] = lower_infill_area_per_combine[0].difference(intersection); // remove thickened area from lower (thickened) layer
                            }
                        }
                    }

                    infill_area_per_combine.push_back(result);
                }
            }
        }
    }
}


}//namespace cura
