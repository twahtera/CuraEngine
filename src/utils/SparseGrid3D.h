/** Copyright (C) 2016 Scott Lenser - Released under terms of the AGPLv3 License */

#ifndef UTILS_SPARSE_GRID_3D_H
#define UTILS_SPARSE_GRID_3D_H

#include "intpoint.h"

#include <cassert>
#include <unordered_map>
#include <vector>

namespace cura {

/*! \brief Sparse grid which can locate spatially nearby elements efficiently in 3D.
 * 
 * \note This is an abstract template class which doesn't have any functions to insert elements.
 * \see SparsePointGrid
 *
 * \tparam ElemT The element type to store.
 */
template<class ElemT>
class SparseGrid3D
{
public:
    using Elem = ElemT;

    /*! \brief Constructs a sparse grid with the specified cell size.
     *
     * \param[in] cell_size The size to use for a cell (square) in the grid.
     *    Typical values would be around 0.5-2x of expected query radius.
     * \param[in] elem_reserve Number of elements to research space for.
     * \param[in] max_load_factor Maximum average load factor before rehashing.
     */
    SparseGrid3D(coord_t cell_size, size_t elem_reserve=0U, float max_load_factor=1.0f);

    static const std::function<bool(const Elem&)> no_precondition;

    /*!
     * Find the nearest element to a given \p query_pt within \p radius.
     *
     * \param[in] query_pt The point for which to find the nearest object.
     * \param[in] radius The search radius.
     * \param[out] elem_nearest the nearest element. Only valid if function returns true.
     * \param[in] precondition A precondition which must return true for an element
     *    to be considered for output
     * \return True if and only if an object has been found within the radius.
     */
    bool getNearest(const Point3 &query_pt, coord_t radius, Elem &elem_nearest,
                    const std::function<bool(const Elem& elem)> precondition = no_precondition) const;

    /*! \brief Process elements from cells that might contain sought after points.
     *
     * Processes elements from cell that might have elements within \p
     * radius of \p query_pt.  Processes all elements that are within
     * radius of query_pt.  May process elements that are up to radius +
     * cell_size from query_pt.
     *
     * \param[in] query_pt The point to search around.
     * \param[in] radius The search radius.
     * \param[in] process_func Processes each element.  process_func(elem) is
     *    called for each element in the cell.
     */
    void processNearby(const Point3 &query_pt, coord_t radius,
                       std::function<bool (const Elem&)>& process_func) const;

    void insert(const Point3 location, const Elem& elem);
    void insert(const std::pair<Point3, Point3> line, const Elem& elem);

    coord_t getCellSize() const;

protected:
    using GridPoint = Point3;
    using grid_coord_t = coord_t;
    using GridMap = std::unordered_multimap<GridPoint, Elem>;

    /*! \brief Process elements from the cell indicated by \p grid_pt.
     *
     * \param[in] grid_pt The grid coordinates of the cell.
     * \param[in] process_func Processes each element.  process_func(elem) is
     *    called for each element in the cell.
     * \return whether to continue processing the next cell
     */
    bool processFromCell(const GridPoint &grid_pt,
                         std::function<bool (const Elem&)>& process_func) const;

    /*!
     * Process each cell along a line segment.
     * 
     * \param line The start and end point of the line segment
     * \param process_func a function executing on each cell along the line; returns true if further cells need to be processed
     */
    void processLine(std::pair<Point3, Point3> line, std::function<bool (const GridPoint&)>& process_func);

    /*! \brief Compute the grid coordinates of a point.
     *
     * \param[in] point The actual location.
     * \return The grid coordinates that correspond to \p point.
     */
    GridPoint toGridPoint(const Point3& point) const;

    /*! \brief Compute the grid coordinate of a print space coordinate.
     *
     * \param[in] coord The actual location.
     * \return The grid coordinate that corresponds to \p coord.
     */
    grid_coord_t toGridCoord(const coord_t& coord) const;

    /*! \brief Compute the lowest point in a grid cell.
     * The lowest point is the point in the grid cell closest to the origin.
     *
     * \param[in] location The grid location.
     * \return The print space coordinates that correspond to \p location.
     */
    Point3 toLowerCorner(const GridPoint& location) const; 

    /*! \brief Compute the lowest coord in a grid cell.
     * The lowest point is the point in the grid cell closest to the origin.
     *
     * \param[in] grid_coord The grid coordinate.
     * \return The print space coordinate that corresponds to \p grid_coord.
     */
    coord_t toLowerCoord(const grid_coord_t& grid_coord) const; 

    /*! \brief Map from grid locations (GridPoint) to elements (Elem). */
    GridMap m_grid;
    /*! \brief The cell (square) size. */
    coord_t m_cell_size;
};



#define SGI_TEMPLATE template<class ElemT>
#define SGI_THIS SparseGrid3D<ElemT>

SGI_TEMPLATE
SGI_THIS::SparseGrid3D(coord_t cell_size, size_t elem_reserve, float max_load_factor)
{
    assert(cell_size > 0U);

    m_cell_size = cell_size;

    // Must be before the reserve call.
    m_grid.max_load_factor(max_load_factor);
    if (elem_reserve != 0U) {
        m_grid.reserve(elem_reserve);
    }
}

SGI_TEMPLATE
typename SGI_THIS::GridPoint SGI_THIS::toGridPoint(const Point3& point)  const
{
    return GridPoint(toGridCoord(point.x), toGridCoord(point.y), toGridCoord(point.z));
}

SGI_TEMPLATE
typename SGI_THIS::grid_coord_t SGI_THIS::toGridCoord(const coord_t& coord)  const
{
    // This mapping via truncation results in the cells with
    // GridPoint.x==0 being twice as large and similarly for
    // GridPoint.y==0.  This doesn't cause any incorrect behavior,
    // just changes the running time slightly.  The change in running
    // time from this is probably not worth doing a proper floor
    // operation.
    return coord / m_cell_size;
}

SGI_TEMPLATE
typename cura::Point3 SGI_THIS::toLowerCorner(const GridPoint& location)  const
{
    return cura::Point3(toLowerCoord(location.x), toLowerCoord(location.y), toLowerCoord(location.z));
}

SGI_TEMPLATE
typename cura::coord_t SGI_THIS::toLowerCoord(const grid_coord_t& grid_coord)  const
{
    // This mapping via truncation results in the cells with
    // GridPoint.x==0 being twice as large and similarly for
    // GridPoint.y==0.  This doesn't cause any incorrect behavior,
    // just changes the running time slightly.  The change in running
    // time from this is probably not worth doing a proper floor
    // operation.
    return grid_coord * m_cell_size;
}

SGI_TEMPLATE
void SGI_THIS::insert(
    const Point3 location,
    const Elem& element)
{
    m_grid.emplace(toGridPoint(location), element);
}

SGI_TEMPLATE
void SGI_THIS::insert(
    const std::pair<Point3, Point3> line,
    const Elem& element)
{
    std::function<bool (const GridPoint&)> process_func = [&element, this](const GridPoint& grid_loc)
        {
            m_grid.emplace(grid_loc, element);
            return true;
        };
    processLine(line, process_func);
}

SGI_TEMPLATE
bool SGI_THIS::processFromCell(
    const GridPoint& grid_pt,
    std::function<bool (const Elem&)>& process_func) const
{
    auto grid_range = m_grid.equal_range(grid_pt);
    for (auto iter = grid_range.first; iter != grid_range.second; ++iter)
    {
        bool continue_ = process_func(iter->second);
        if (!continue_)
        {
            return false;
        }
    }
    return true;
}

SGI_TEMPLATE
void SGI_THIS::processNearby(const Point3& query_pt, coord_t radius,
                             std::function<bool (const Elem&)>& process_func) const
{
    Point3 min_loc = query_pt - Point3(radius, radius, radius);
    Point3 max_loc = query_pt + Point3(radius, radius, radius);

    GridPoint min_grid = toGridPoint(min_loc);
    GridPoint max_grid = toGridPoint(max_loc);

    for (coord_t grid_y = min_grid.y; grid_y <= max_grid.y; ++grid_y)
    {
        for (coord_t grid_x = min_grid.x; grid_x <= max_grid.x; ++grid_x)
        {
            for (coord_t grid_z = min_grid.z; grid_z <= max_grid.z; ++grid_z)
            {
                GridPoint grid_pt(grid_x, grid_y, grid_z);
                bool continue_ = processFromCell(grid_pt, process_func);
                if (!continue_)
                {
                    return;
                }
            }
        }
    }
}

SGI_TEMPLATE
void SGI_THIS::processLine(const std::pair<Point3, Point3> line, std::function<bool (const GridPoint&)>& process_func)
{
    Point3 diff = line.second - line.first;
    GridPoint start = toGridPoint(line.first);
    GridPoint end = toGridPoint(line.second);

    auto sign = [](int in)
        {
            return (in > 0)? 1 : (in < 0)? -1 : 0;
        };

    auto sign0 = [](int in)
        {
            return (in > 0)? 1 : 0;
        };

    bool continue_ = process_func(start);
    if (!continue_)
    {
        return;
    }
    if (diff == Point3(0,0,0))
    {
        return;
    }
    int x_dir = sign(diff.x);
    int y_dir = sign(diff.y);
    int z_dir = sign(diff.z);
    Point3 last_cell_crossing = line.first;
    for (GridPoint now = start; now != end;)
    {
        float best_progress;
        int best_axis = -1;
        if (x_dir)
        {
            coord_t next_x_coord = toLowerCoord(now.x + sign0(diff.x));
            if (next_x_coord == last_cell_crossing.x)
            {
                next_x_coord = toLowerCoord(now.x + sign(diff.x));
            }
            float next_x_progress = (next_x_coord - line.first.x) / (line.second.x - line.first.x);
            best_progress = next_x_progress;
            best_axis = 0;
        }
        if (y_dir)
        {
            coord_t next_y_coord = toLowerCoord(now.y + sign0(diff.y));
            if (next_y_coord == last_cell_crossing.y)
            {
                next_y_coord = toLowerCoord(now.y + sign(diff.y));
            }
            float next_y_progress = (next_y_coord - line.first.y) / (line.second.y - line.first.y);
            if (next_y_progress < best_progress || best_axis == -1)
            {
                best_progress = next_y_progress;
                best_axis = 1;
            }
        }
        if (z_dir)
        {
            coord_t next_z_coord = toLowerCoord(now.z + sign0(diff.z));
            if (next_z_coord == last_cell_crossing.z)
            {
                next_z_coord = toLowerCoord(now.z + sign(diff.z));
            }
            float next_z_progress = (next_z_coord - line.first.z) / (line.second.z - line.first.z);
            if (next_z_progress < best_progress || best_axis == -1)
            {
                best_progress = next_z_progress;
                best_axis = 2;
            }
        }
        Point3 next_x_crossing = line.first + diff * best_progress;
        GridPoint next = toGridPoint(next_x_crossing);
        assert(next != now && (next - now).vSize2() == 1 && "We must take step of size one!");
        bool continue_ = process_func(now);
        if (!continue_)
        {
            return;
        }
    }
}

SGI_TEMPLATE
const std::function<bool(const typename SGI_THIS::Elem &)>
    SGI_THIS::no_precondition =
    [](const typename SGI_THIS::Elem &)
    {
        return true;
    };

SGI_TEMPLATE
bool SGI_THIS::getNearest(
    const Point3& query_pt, coord_t radius, Elem& elem_nearest,
    const std::function<bool(const Elem& elem)> precondition) const
{
    bool found = false;
    int64_t best_dist2 = static_cast<int64_t>(radius) * radius;
    std::function<bool (const Elem&)> process_func =
        [&query_pt, &elem_nearest, &found, &best_dist2, &precondition](const Elem& elem)
        {
            if (!precondition(elem))
            {
                return true;
            }
            int64_t dist2 = vSize2(elem.point - query_pt);
            if (dist2 < best_dist2)
            {
                found = true;
                elem_nearest = elem;
                best_dist2 = dist2;
            }
            return true;
        };
    processNearby(query_pt, radius, process_func);
    return found;
}

SGI_TEMPLATE
coord_t SGI_THIS::getCellSize() const
{
    return m_cell_size;
}

#undef SGI_TEMPLATE
#undef SGI_THIS

} // namespace cura

#endif // UTILS_SPARSE_GRID_3D_H
