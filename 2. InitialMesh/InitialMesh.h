#ifndef LIBRARY_INITIAL_MESH
#define LIBRARY_INITIAL_MESH

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <string>
#include <iostream>
#include <fstream>

//#include <algorithm>
#include <numbers>
//using std::numbers::pi;

class InitialMesh {

private:

	static constexpr double PI = std::numbers::pi;
	static constexpr double DOUBLEMAX = std::numeric_limits<double>::max();
	
		
	struct Point
	{
		//Point() = default;	
		Point(double xx, double yy, double rr) :x{ xx }, y{ yy }, r{ rr } {};
		double x{ 0.0 }, y{ 0.0 }, r{ 0.0 };
	};
	using Points = std::vector<Point>;	


	struct Triangle
	{
		//Triangle() = default;
		Triangle(std::size_t aa, std::size_t bb, std::size_t cc)
			: a{ aa }, b{ bb }, c{ cc } {}
		std::size_t a{ 0 }, b{ 0 }, c{ 0 };
	};
	using Triangles = std::vector<Triangle>;
	
		
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

	struct EdgeP
	{
		//EdgeP() = default;
		EdgeP(std::size_t NT, std::size_t C)
			:nt{ NT }, c{ C } {}
		std::size_t nt{ 0 }, c{ 0 };
	};
	using EdgePs = std::list<EdgeP>;

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
	
	using EdgeEdgePss = std::unordered_map<Edge, EdgePs, EdgeHash, EdgeEqual>;


	//using Edges = std::unordered_set<Edge, EdgeHash, EdgeEqual>;

	using BoundaryList = std::list<std::size_t>;
	using BoundaryLists = std::vector<std::list<std::size_t>>;
	
	
	double getAngle(double xa, double ya, double xb, double yb, double xi, double yi) const
	{
		return std::acos(((xa - xi) * (xb - xi) + (ya - yi) * (yb - yi)) /
			(std::hypot(xa - xi, ya - yi) * std::hypot(xb - xi, yb - yi)));
	}

	double getHeightOfTriangle(double xa, double ya, double xb, double yb, double xi, double yi) const
	{
		return std::abs((xa - xi) * (yb - yi) - (xb - xi) * (ya - yi))
			/ std::hypot(xb - xa, yb - ya);
	}



public:

	// Data input
	void input();

	// mesh data initialization. 
	void initialization();

	//extent is used to control the initial radii of point bubbles 
	//(the value is suggested to be larger than 2.0) 
	// ratio: ratio of neighboring element sizes, the value should be larger than 1.0.
	//ratio_boundary: control the distance of neighbouring boundaries,
	//the larger value refelct larger distance of neighbouring boundaries,
	// the value is suggested to be larger than 3.0.
	void meshGeneration(double extent, double ratio_Distence, double ratio_Radius, double ratio_Gradient);


	//mesh output 
	void output() const;

private:


	Points ps;	Triangles ts; EdgeEdgePss eepss;
	BoundaryLists bls;

	
	// erase Triange Edges and EdgePs from EdgeEdgePss
	void eraseTriangeEdges(std::size_t i);

	// insert Triange Edges and EdgePs to EdgeEdgePss
	void insertTriangeEdges(std::size_t i);

	//extent is used to control the initial radii of point bubbles 
	//(the value is suggested to be larger than 2.0) 
	void setRadiiOfPointBubblesAccordingToRefinementExtent(double extent);
	
	void updateRadiiOfPointBubblesAccordingToDistanceRatio(double ratio_Distance);

	
    void updateRadiiOfPointBubblesAccordingToGradientRatio(double ratio);

	double getPredefinedMinimumRadiusAccordingToRadiusRatio(double ratio_Radius) const;

	void refineBoundaryEdgesAccordingToGradientRatio(double ratio_Gradient);	

	
	using ThetaEs = std::multimap<double, const Edge>;
	using EIts = std::unordered_map<Edge, ThetaEs::iterator, EdgeHash, EdgeEqual>;
	void insert(ThetaEs& thetaes, EIts& eits, const Edge& e);
	void erase(ThetaEs& thetaes, EIts& eits, const Edge& e);
	void updateConnectionsOfPoints();


	void updateRadiiOfPointBubblesAccordingToDistanceRatioAndMinimumRadius(double ratio_Distence, double rmin);

	
	void adjuestBoundaryPointsPositions();


	
};

#endif



//using ETHETAS = std::unordered_map<Edge, double, EdgeHash, EdgeEqual>;

//void insert(ETHETAS& ethetas, const Edge& e);