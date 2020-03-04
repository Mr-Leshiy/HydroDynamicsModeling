#include <triangulation.h>

#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(triangulation_tests)

BOOST_AUTO_TEST_CASE(point_struct_test)
{
    TPoint point1(4, 4, 4);
    TPoint point2(2, 2, 2);

    point1 += point2;

    TPoint expected(6, 6, 6);

    BOOST_CHECK(point1.x() == expected.x() && point1.y() == expected.y() && point1.z() == expected.z());

    point1 -= point2;

    expected = TPoint(4, 4, 4);

    BOOST_CHECK(point1.x() == expected.x() && point1.y() == expected.y() && point1.z() == expected.z());

    point1 /= 2;

    expected = TPoint(2, 2, 2);

    BOOST_CHECK(point1.x() == expected.x() && point1.y() == expected.y() && point1.z() == expected.z());

    point1 *= 4;

    expected = TPoint(8, 8, 8);

    BOOST_CHECK(point1.x() == expected.x() && point1.y() == expected.y() && point1.z() == expected.z());
}

BOOST_AUTO_TEST_CASE(simple_test)
{
    P3DT3::Point_3 point1(1, 2, 3);
    P3DT3::Point_3 point2(2, 3, 4);
    P3DT3::Point_3 point3(1, 2, 2);
    P3DT3::Point_3 point4(3, 2, 3);
    P3DT3::Point_3 point5(-3, 2, 3);

    BOOST_CHECK(point1 == point2);

    P3DT3::Tetrahedron tet1(point1, point2, point3, point4);
    P3DT3::Tetrahedron tet2(point1, point2, point3, point4);

    BOOST_CHECK(tet1[0].x() == 1);

    P3DT3::Tetrahedron& tet3 = tet1;
    P3DT3::Tetrahedron& tet4 = tet2;

    BOOST_CHECK(tet3 == tet4);
}

BOOST_AUTO_TEST_CASE(simple_test2)
{
    double box_size = 10;

    P3DT3::Point_3 point1(-0.02, 0.085, 0.308); 
    P3DT3::Point_3 point2(0.358, -0.325, 0.246);
    P3DT3::Point_3 point3(-0.135, -0.408, -0.485);
    P3DT3::Point_3 point4(-0.495, -0.380, -0.054);
    P3DT3::Point_3 point5(0.107, 0.101, 0.071);

    std::vector<P3DT3::Point_3> points;

    P3DT3 delaunay_triangulation;


    delaunay_triangulation.insert(point1);
    delaunay_triangulation.insert(point2);
    delaunay_triangulation.insert(point3);
    delaunay_triangulation.insert(point4);
    delaunay_triangulation.insert(point5);

    point1 += Gt::Point_3(1, 2, 3);

}


BOOST_AUTO_TEST_SUITE_END()