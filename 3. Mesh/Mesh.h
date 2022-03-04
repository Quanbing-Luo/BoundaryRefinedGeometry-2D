#ifndef LIBRARY_MESH
#define LIBRARY_MESH

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <numbers>
#include <string>
#include <fstream>
#include <algorithm>


class Mesh
{
private:
	static constexpr double PI = std::numbers::pi;
	//using Numbers = std::vector<std::size_t>;
	
	struct Point
	{
		//Point() = default;	
		Point(double xx, double yy, double rr) :x{ xx }, y{ yy }, r{rr} {};
		double x{ 0.0 }, y{ 0.0 }, r{ 0.0 };
	};
	using Points = std::vector<Point>;

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


	struct Triangle
	{
		//Triangle() = default;
		Triangle(std::size_t aa, std::size_t bb, std::size_t cc)
			: a{ aa }, b{ bb }, c{ cc }, ec3{std::make_pair(Edge(a,b),c),
			std::make_pair(Edge(b,c),a) ,std::make_pair(Edge(c,a),b) } { }

		std::size_t a{ 0 }, b{ 0 }, c{ 0 };		

		std::array<std::pair<Edge, std::size_t>, 3> ec3;

		void setParameters(const Points& ps);

		// circumcenter and radius of Circumcircle
		double xcc{ 0.0 }, ycc{ 0.0 }, rcc{ 0.0 };

		// the center of gravity and length scale (radius)
		double x{ 0.0 }, y{ 0.0 }, r{ 0.0 };

		//area
		double area{ 0.0 }, area_bubble{ 0.0 };

		//mesh density
		double density{ 0.0 };
	};

	using Triangles = std::vector<Triangle>;



	struct EdgeP
	{
		EdgeP() = default;
		
		// number of edge 
		std::size_t sz{ 0 };
		// position of triangle and point c
		std::size_t n1{ 0 }, c1{ 0 };
		std::size_t n2{ 0 }, c2{ 0 };	

		double length{ 0.0 };
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

	using EdgeEdgePs = std::unordered_map<Edge, EdgeP, EdgeHash, EdgeEqual>;
	
public:
	
	// data input
	void input();

	// data initialization. 
	void initialization();

	// delaunay triangulation. 
	void meshGeneration();

	//mesh output 
	void output() const;


private:

	Points ps;	Triangles ts;  EdgeEdgePs eeps;


	using Numbers = std::vector<std::size_t>; 
	using Edges = std::unordered_set<Edge,  EdgeHash, EdgeEqual>;
	using CavityPair = std::pair<Numbers, Edges>;
	CavityPair getCavity(std::size_t nt) const;

	void insertTriange(std::size_t n);

	void eraseTriange(std::size_t n);

};

#endif
