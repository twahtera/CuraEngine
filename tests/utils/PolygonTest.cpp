//Copyright (c) 2015 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include "PolygonTest.h"

#include "../src/utils/SVG.h" // TODO temp

namespace cura
{
    CPPUNIT_TEST_SUITE_REGISTRATION(PolygonTest);

void PolygonTest::setUp()
{
    test_square.emplace_back(0, 0);
    test_square.emplace_back(100, 0);
    test_square.emplace_back(100, 100);
    test_square.emplace_back(0, 100);


    pointy_square.emplace_back(0, 0);
    pointy_square.emplace_back(47, 0);
    pointy_square.emplace_back(50, 80);
    pointy_square.emplace_back(53, 0);
    pointy_square.emplace_back(100, 0);
    pointy_square.emplace_back(100, 100);
    pointy_square.emplace_back(55, 100);
    pointy_square.emplace_back(50, 180);
    pointy_square.emplace_back(45, 100);
    pointy_square.emplace_back(0, 100);

    pointy_square2.emplace_back(0, 0);
    pointy_square2.emplace_back(400, 0);
    pointy_square2.emplace_back(470, 800);
    pointy_square2.emplace_back(530, 800);
    pointy_square2.emplace_back(600, 0);
    pointy_square2.emplace_back(1000, 0);
    pointy_square2.emplace_back(1000, 1000);
    pointy_square2.emplace_back(550, 1000);
    pointy_square2.emplace_back(500, 1800);
    pointy_square2.emplace_back(450, 1000);
    pointy_square2.emplace_back(0, 1000);

    triangle.emplace_back(100, 0);
    triangle.emplace_back(300, 0);
    triangle.emplace_back(200, 100);

    clipper_bug.emplace_back(107347, 120836);
    clipper_bug.emplace_back(107309, 120910);
    clipper_bug.emplace_back(107158, 120960);
    clipper_bug.emplace_back(106760, 120839);
    clipper_bug.emplace_back(106570, 120831);

    PolygonRef outline = complex.newPoly();
    outline.emplace_back(0, 0);
    outline.emplace_back(400, 0);
    outline.emplace_back(470, 300);
    outline.emplace_back(530, 300);
    outline.emplace_back(600, 0);
    outline.emplace_back(1000, 0);
    outline.emplace_back(1000, 1000);
    outline.emplace_back(1300, 1050);
    outline.emplace_back(1300, 1150);
    outline.emplace_back(1000, 1200);
    outline.emplace_back(1000, 1500);
    outline.emplace_back(550, 1500);
    outline.emplace_back(500, 1800);
    outline.emplace_back(450, 1500);
    outline.emplace_back(0, 1500);
    outline.emplace_back(0, 1200);
    outline.emplace_back(300, 1100);
    outline.emplace_back(0, 1000);
    PolygonRef hole = complex.newPoly();
    hole.emplace_back(305, 500);
    hole.emplace_back(300, 905);
    hole.emplace_back(705, 900);
    hole.emplace_back(700, 500);
    
}

void PolygonTest::tearDown()
{
    //Do nothing.
}


void PolygonTest::polygonOffsetTest()
{
    Polygons test_squares;
    test_squares.add(test_square);
    Polygons expanded = test_squares.offset(25);
    int64_t expanded_length = expanded.polygonLength();

    Polygons square_hole;
    PolygonRef square_inverted = square_hole.newPoly();
    for (int i = test_square.size() - 1; i >= 0; i--)
    {
        square_inverted.add(test_square[i]);
    }
    Polygons contracted = square_hole.offset(25);
    int64_t contracted_length = contracted.polygonLength();

    CPPUNIT_ASSERT_MESSAGE("Offset on outside poly is different from offset on inverted poly!", std::abs(expanded_length - contracted_length) < 5);
}

void PolygonTest::polygonOffsetBugTest()
{
    Polygons ins = complex;
    Polygons in_between = ins.cornerBasedOffset(-100, 100, -10);
    in_between.removeSmallAreas(1);
    Polygons in_between2 = in_between.cornerBasedOffset(-100, 100, -10);
    Polygons result = in_between2.cornerBasedOffset(-100, 100, -10);
    {
        AABB aabb(ins);
        for (PolygonRef poly : result)
            for (const Point p : poly)
                aabb.include(p);
        SVG svg("poly_test.html", aabb, Point(1024 * 1, 1024 * 1));
        svg.writePolygons(ins, SVG::Color::BLACK);
        svg.writePoints(ins, true);
        svg.writePolygons(result, SVG::Color::BLUE);
    }
    std::cerr << "done\n";
}


void PolygonTest::cornerBasedOffsetTest()
{
    Polygons polys;
    polys.add(clipper_bug);
    Polygons offsetted = polys.offset(-20);

    for (PolygonRef poly : offsetted)
    {
        for (Point& p : poly)
        {
            CPPUNIT_ASSERT_MESSAGE("Polygon offset moved point the wrong way!", polys.inside(p));
        }
    }
}


void PolygonTest::isOutsideTest()
{
    Polygons test_triangle;
    test_triangle.add(triangle);

    CPPUNIT_ASSERT_MESSAGE("Left point is calculated as inside while it's outside!", !test_triangle.inside(Point(0, 100)));
    CPPUNIT_ASSERT_MESSAGE("Middle left point is calculated as inside while it's outside!", !test_triangle.inside(Point(100, 100)));
    CPPUNIT_ASSERT_MESSAGE("Middle right point is calculated as inside while it's outside!", !test_triangle.inside(Point(300, 100)));
    CPPUNIT_ASSERT_MESSAGE("Right point is calculated as inside while it's outside!", !test_triangle.inside(Point(500, 100)));
    CPPUNIT_ASSERT_MESSAGE("Above point is calculated as inside while it's outside!", !test_triangle.inside(Point(100, 200)));
    CPPUNIT_ASSERT_MESSAGE("Below point is calculated as inside while it's outside!", !test_triangle.inside(Point(100, -100)));
}

void PolygonTest::isInsideTest()
{
    Polygons test_polys;
    PolygonRef poly = test_polys.newPoly();
    poly.add(Point(82124,98235));
    poly.add(Point(83179,98691));
    poly.add(Point(83434,98950));
    poly.add(Point(82751,99026));
    poly.add(Point(82528,99019));
    poly.add(Point(81605,98854));
    poly.add(Point(80401,98686));
    poly.add(Point(79191,98595));
    poly.add(Point(78191,98441));
    poly.add(Point(78998,98299));
    poly.add(Point(79747,98179));
    poly.add(Point(80960,98095));

    CPPUNIT_ASSERT_MESSAGE("Inside point is calculated as being outside!", test_polys.inside(Point(78315, 98440)));

}




}