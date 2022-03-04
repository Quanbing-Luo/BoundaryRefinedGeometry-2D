#ifndef LIBRARY_OptimizedMesh
#define LIBRARY_OptimizedMesh
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <utility>
//#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <numbers>

class OptimizedMesh
{
private:

	static constexpr double PI = std::numbers::pi;

	struct Edge
	{
		//Edge() = default;
		Edge(std::size_t aa, std::size_t bb)
		{
			if (aa < bb) { a = aa; b = bb; }
			else { a = bb; b = aa; }
		}
		std::size_t a{ 0 }, b{ 0 };
	};

	struct EdgeHash
	{
		std::size_t operator()(const Edge& e) const
		{
			return std::hash<std::size_t>{}(e.a) ^ std::hash<std::size_t>{}(e.b);
		}
	};

	struct EdgeEqual {
		bool operator()(const Edge& e1, const Edge& e2) const
		{
			return (e1.a == e2.a) && (e1.b == e2.b);
		}
	};


	struct EdgeP
	{
		//EdgeP() = default;
		EdgeP(std::size_t NT, std::size_t C)
			:nt{ NT }, c{ C } {}
		std::size_t nt{ 0 }, c{ 0 };
	};

	using Edges = std::unordered_set<Edge, EdgeHash, EdgeEqual>;
	using EdgePs = std::list<EdgeP>;
	using EdgeEdgePss = std::unordered_map<Edge, EdgePs, EdgeHash, EdgeEqual>;


	// Point 
	struct Point
	{
		Point() = default;

		Point(double xx, double yy, double rr)
			:x{ xx }, y{ yy }, r{ rr } {}

		double x{ 0.0 }, y{ 0.0 }, r{ 0.0 };
		double dx{ 0.0 }, dy{ 0.0 }, dr{ 0.0 };
		bool isatboundary{ false };

		//Around Points
		std::unordered_set<size_t> aaps;
	};
	using Points = std::vector<Point>;

	//Triangle and Triangles 
	struct Triangle
	{
		//Triangle() = default;
		Triangle(std::size_t ii, std::size_t aa, std::size_t bb, std::size_t cc)
			: a{ aa }, b{ bb }, c{ cc },
			eep3{ std::make_pair(Edge(a,b),EdgeP(ii,c)),
				 std::make_pair(Edge(b,c),EdgeP(ii,a)),
				 std::make_pair(Edge(c,a),EdgeP(ii,b)) } { }	

		std::size_t a{ 0 }, b{ 0 }, c{ 0 };

		std::array<std::pair<Edge, EdgeP>, 3> eep3;
	};
	using Triangles = std::vector<Triangle>;


public:

	// mesh data input 
	void input();

	//initialization
	void initialization();

	//optimization
	void optimization();

	//mesh output 
	void output() const;
	

private:

	Points ps;	Triangles ts;   EdgeEdgePss eepss;


	// erase Triange Edges and EdgePs from EdgeEdgePss
	void eraseTriangeEdges(std::size_t n);

	// insert Triange Edges and EdgePs to EdgeEdgePss
	void insertTriangeEdges(std::size_t n);


	//radii calculation
	void radii_calculation();

	//new positions calculation
	void new_positions_calculation();

	double getRelativeLength(std::size_t a, std::size_t b, std::size_t c1, std::size_t c2) const;
	void updatePointsConnections();

};

#endif