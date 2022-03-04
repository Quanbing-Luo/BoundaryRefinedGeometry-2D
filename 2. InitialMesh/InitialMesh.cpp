#include "InitialMesh.h"

// Data input
void InitialMesh::input()
{
	std::ifstream ist;
	ist.open("Geometry.txt");
	std::string ss;	

	// input Points
	std::getline(ist, ss);
	std::size_t np{ 0 }; ist >> np;
	for (std::size_t i = 0; i < np; i++)
	{
		static std::size_t n{ 0 };
		static double x{ 0.0 }, y{ 0.0 };

		ist >> n >> x >> y ;
		ps.emplace_back(x, y, DOUBLEMAX);
	}
	ist >> std::ws;
	
	// input Triangles 
	std::getline(ist, ss);
	std::size_t nt{ 0 }; ist >> nt;
	for (std::size_t i = 0; i < nt; i++)
	{
		static std::size_t n{ 0 };
		static std::size_t a{ 0 }, b{ 0 }, c{ 0 };

		ist >> n >> a >> b >> c ;
		ts.emplace_back(a, b, c);
	}
	ist >> std::ws;	

	ist.close();
}


void InitialMesh::output() const
{
	//const std::string name{ };
	std::ofstream ost;
	ost.open("InitialMesh.txt");

	// output Points 
	ost << "Points\n";
	ost << ps.size() << '\n';
	for (std::size_t i = 1; auto & p : ps)
	{
		ost << std::scientific << i << '\t'
			<< p.x << '\t' << p.y << '\t' << p.r  << '\n';
		i++;
	}

	// output Triangles 
	ost << "Triangles\n";
	ost << ts.size() << '\n';
	for (std::size_t i = 1; auto & t : ts)
	{
		ost << std::scientific << i << '\t'
			<< t.a << '\t' << t.b << '\t' << t.c << '\n';
		i++;
	}

	ost.close();
}


void InitialMesh::insertTriangeEdges(std::size_t i)
{
	auto& [a, b, c] { ts[i]};
	eepss[Edge(a, b)].emplace_back(EdgeP(i, c));
	eepss[Edge(b, c)].emplace_back(EdgeP(i, a));
	eepss[Edge(c, a)].emplace_back(EdgeP(i, b));
}


void InitialMesh::eraseTriangeEdges(std::size_t i)
{
	auto& [a, b, c] { ts[i] };

	std::vector<Edge> e3{ Edge(a,b), Edge(a,c),Edge(b,c) };
	std::vector<EdgeP> ep3{ EdgeP(i,c), EdgeP(i,b),EdgeP(i,a) };

	for (std::size_t i = 0; i < 3; i++)
	{
		auto& e{ e3[i] };
		auto& ep{ ep3[i] };

		auto iteeps = eepss.find(e);
		auto& eps{ iteeps->second };
		for (auto it = eps.begin(); it != eps.end(); it++)
		{
			if (it->nt == ep.nt && it->c == ep.c)
			{
				eps.erase(it);
				break;
			}
		}

		if (eps.size() == 0)
			eepss.erase(iteeps);
	}
}



void InitialMesh::initialization()
{
	for (std::size_t i=0; auto& t : ts)
	{		
		insertTriangeEdges(i);
		i++;
	}

	for (auto& [e, eps] : eepss)
	{
		if (eps.size() == 1)			
			bls.emplace_back(BoundaryList{ e.a, e.b });
	}				
}


//
void InitialMesh::setRadiiOfPointBubblesAccordingToRefinementExtent(double extent)
{
	// set radii with extent 

	for (const auto& [e, eps] : eepss)
	{
		auto& [a, b] {e};
		auto& [xa, ya, ra] {ps[a]};
		auto& [xb, yb, rb] {ps[b]};
		double rr{ std::hypot(xb - xa, yb - ya) / extent };
		if (rr < ra) ra = rr;
		if (rr < rb) rb = rr;		
	}
}


void InitialMesh::updateRadiiOfPointBubblesAccordingToDistanceRatio(double ratio_Distance)
{
	for (const auto& [e, eps] : eepss)
	{
		if (eps.size() == 1)
		{
			auto& [a, b] {e};
			std::size_t c{ eps.front().c };

			auto& [xa, ya, ra] {ps[a]};
			auto& [xb, yb, rb] {ps[b]};
			auto& [xc, yc, rc] {ps[c]};

			double theta_a = getAngle(xc, yc, xb, yb, xa, ya);
			double theta_b = getAngle(xc, yc, xa, ya, xb, yb);

			if (theta_a < (PI / 2.0 + 1e-8) && theta_b < (PI / 2.0 + 1e-8))
			{
				double rtemp = getHeightOfTriangle(xa, ya, xb, yb, xc, yc) / ratio_Distance;
				if (rtemp < rc)
					rc = rtemp;
			}
				
		}
	}
}


void InitialMesh::updateRadiiOfPointBubblesAccordingToGradientRatio(double ratio_Gradient)
{	
	std::vector<std::unordered_set<std::size_t>> nbss(ps.size());
	for (const auto& [e, eps] : eepss)
	{
		auto& [a, b] {e};
		nbss[a].insert(b);
		nbss[b].insert(a);
	}

	using RIs = std::multimap<double, std::size_t>;
	using IIts = std::unordered_map<std::size_t, RIs::iterator>;
	RIs ris; IIts iits;

	for (std::size_t i = 0; auto & p:ps)
	{
		auto it = ris.emplace(p.r, i);
		iits.emplace(i, it);
		i++;
	}

	while (ris.size() > 0)
	{
		auto itmin = ris.begin();
		auto a{ itmin->second };

		//delete min
		ris.erase(itmin); 	iits.erase(a);

		auto& [xa, ya, ra] { ps[a] };
		for (auto& b : nbss[a])
		{
			auto& [xb, yb, rb] { ps[b] };
			if (rb > ra)
			{
				double rtemp{ ra + (ratio_Gradient - 1.0) / (ratio_Gradient + 1.0) * std::hypot(xb - xa, yb - ya) };
				if (rb > rtemp)
				{
					rb = rtemp;

					//update radius b
					auto& itold = iits.at(b);  		ris.erase(itold);
					auto itnew = ris.emplace(rb, b);  itold = itnew;
				}
			}
		}
	}
}




double InitialMesh::getPredefinedMinimumRadiusAccordingToRadiusRatio(double ratio_Radius) const
{
	double rmin = DOUBLEMAX;
	for (const auto& [x, y, r] : ps)
	{
		if (r < rmin)
			rmin = r;
	}

	return rmin / ratio_Radius;
}



void InitialMesh::refineBoundaryEdgesAccordingToGradientRatio(double ratio_Gradient)
{	
	for (auto& bl : bls)
	{
		double sumradii{ 0.0 };
		std::multimap<double, BoundaryList::iterator> ratioits;
		for (auto ita = bl.begin(), itb = std::next(ita); itb != bl.end(); ita++, itb++)
		{
			auto& [xa, ya, ra] {ps[*ita]}; auto& [xb, yb, rb] {ps[*itb]};

			auto ratio = (ra + rb)/ std::hypot(xb - xa, yb - ya);
			ratioits.emplace(ratio, ita);
			sumradii += ra + rb;
		}

		std::size_t aa{ bl.front() }, bb{ bl.back() };
		auto& [xaa, yaa, raa] {ps[aa]};	auto& [xbb, ybb, rbb] {ps[bb]};
		double length = std::hypot(xbb - xaa, ybb - yaa);

		while (sumradii < length)
		{
			auto itratiomin = ratioits.begin();
			auto ita{ itratiomin->second }; auto  itb{ std::next(ita) };
			std::size_t a{ *ita }, b{ *itb };
			auto emin{ Edge(a,b) };
			auto& epsmin = eepss.at(emin);
			auto [nt, c] {epsmin.front()};

			auto [xa, ya, ra] { ps[a] }; auto [xb, yb, rb] { ps[b] };
			double xi{ (xa + xb) / 2.0 }, yi{ (ya + yb) / 2.0 };

			double rtemp1{ ra + (ratio_Gradient - 1.0) / (ratio_Gradient + 1.0) * std::hypot(xi - xa, yi - ya) };
			double rtemp2{ rb + (ratio_Gradient - 1.0) / (ratio_Gradient + 1.0) * std::hypot(xi - xb, yi - yb) };
			double ri{ std::min(rtemp1,rtemp2) };


			double abs1{ std::abs(length - sumradii) };
			double abs2{ std::abs(length - sumradii - 2.0 * ri) };
			if (abs1 > abs2)
			{
				eraseTriangeEdges(nt);

				ratioits.erase(itratiomin);


				ps.emplace_back(xi, yi, ri);
				std::size_t i{ ps.size() - 1 };

				ts[nt] = Triangle( a, i, c);
				std::size_t newnt = ts.size();
				ts.emplace_back(i, b, c);

				insertTriangeEdges(nt);
				insertTriangeEdges(newnt);

				auto  iti = bl.insert(itb, i);

				double ratioa = (ra + ri) / std::hypot(xi - xa, yi - ya) ;
				ratioits.emplace(ratioa, ita);
				double ratioi = (rb + ri) / std::hypot(xb - xi, yb - yi) ;
				ratioits.emplace(ratioi, iti);

				sumradii += 2.0 * ri;
				//output();
				//continue;
			}
			else
				break;
		}

		//output();
	}	
}

void InitialMesh::insert(ThetaEs& thetaes, EIts& eits, const Edge& e)
{
	auto& eps = eepss.at(e);
	if (eps.size() == 2)
	{
		auto& [a, b] {e};
		std::size_t c1 = eps.front().c;
		std::size_t c2 = eps.back().c;
		auto& [xa, ya, ra] {ps[a]}; auto& [xb, yb, rb] {ps[b]};
		auto& [xc1, yc1, rc1] {ps[c1]}; auto& [xc2, yc2, rc2] {ps[c2]};

		double theta = 2.0 * PI - getAngle(xa, ya, xb, yb, xc1, yc1) - getAngle(xa, ya, xb, yb, xc2, yc2);;
		if (theta < PI - 1e-8)
		{
			auto it = thetaes.emplace(theta, e);
			eits.emplace(e, it);
		}
	}
}

void InitialMesh::erase(ThetaEs& thetaes, EIts& eits, const Edge& e)
{
	auto it = eits.find(e);
	if (it != eits.end())
	{
		thetaes.erase(it->second);
		eits.erase(it);
	}
}

void InitialMesh::updateConnectionsOfPoints()
{
	ThetaEs thetaes;  EIts eits;

	for (auto& [e, eps] : eepss)
		insert(thetaes, eits, e);

	while (thetaes.size() > 0)
	{
		auto itmin = thetaes.begin();
		auto& emin{ itmin->second };	auto& epsmin = eepss.at(emin);
		auto [a, b] {emin};
		auto [nt1, c1] { epsmin.front() };
		auto [nt2, c2] { epsmin.back()};


		//delete olds 
		eraseTriangeEdges(nt1);
		eraseTriangeEdges(nt2);

		erase(thetaes, eits, emin);
		erase(thetaes, eits, Edge(a, c1));
		erase(thetaes, eits, Edge(b, c1));
		erase(thetaes, eits, Edge(a, c2));
		erase(thetaes, eits, Edge(b, c2));


		//add news
		ts[nt1] = Triangle(c1, c2, a);
		ts[nt2] = Triangle(c1, c2, b);
	
		insertTriangeEdges(nt1);
		insertTriangeEdges(nt2);

		insert(thetaes, eits, Edge(a, c1));
		insert(thetaes, eits, Edge(a, c2));
		insert(thetaes, eits, Edge(b, c1));
		insert(thetaes, eits, Edge(b, c2));
		//insert(ethetas, Edge(c1, c2));	

		//output();
	}	
}


void InitialMesh::updateRadiiOfPointBubblesAccordingToDistanceRatioAndMinimumRadius(double ratio_Distance, double rmin)
{
	for (const auto& [e, eps] : eepss)
	{
		if (eps.size() == 1)
		{
			auto& [a, b] {e};
			std::size_t c{ eps.front().c };

			if (eepss.at(Edge(a, c)).size() == 1 || eepss.at(Edge(b, c)).size() == 1)
				continue;
			else
			{
				auto& [xa, ya, ra] {ps[a]};
				auto& [xb, yb, rb] {ps[b]};
				auto& [xc, yc, rc] {ps[c]};

				double theta_a = getAngle(xc, yc, xb, yb, xa, ya);
				double theta_b = getAngle(xc, yc, xa, ya, xb, yb);

				if (theta_a < (PI / 2.0 + 1e-8) && theta_b < (PI / 2.0 + 1e-8))
				{
					double rtemp = getHeightOfTriangle(xa, ya, xb, yb, xc, yc) / ratio_Distance;
					if(rtemp<rmin)
						rc = rmin;
					else if (rtemp < rc)
						rc = rtemp;
				}
			}
		}
	}
}





void InitialMesh::adjuestBoundaryPointsPositions()
{	
	for (auto& bl : bls)
	{
		if (bl.size() == 2)
			continue;
		
		std::size_t a{ bl.front() };
		std::size_t b{ bl.back() };
		const auto& [xa, ya, ra] {ps[a]};
		const auto& [xb, yb, rb] {ps[b]};
		double ERR = std::hypot(xb - xa, yb - ya) / ((bl.size() - 1) * 1e4);

		std::vector<double> dxs(bl.size() - 2);
		std::vector<double> dys(bl.size() - 2);			

		for (std::size_t n = 0; n < 1e10; n++)
		{
			double err = 0.0;
			auto itbegin = std::next(bl.begin()), itend = std::prev(bl.end());
			
			std::size_t i = 0;
			for (auto itc = itbegin; itc != itend; itc++)
			{
				auto itl = std::prev(itc), itr = std::next(itc);				
				std::size_t l{ *itl }, c{ *itc }, r{ *itr };				
				auto& [xl, yl, rl] { ps[l] }; 
				auto& [xc, yc, rc] { ps[c] }; 
				auto& [xr, yr, rr] { ps[r] };

				double lengthl = std::hypot(xl - xc, yl - yc);
				double lengthr = std::hypot(xr - xc, yr - yc);
				double rrl = rl + rc; double rrr = rr + rc;
				double fl = lengthl - rrl; double fr = lengthr - rrr;

				double nxl = (xl - xc) / lengthl; double nxr = (xr - xc) / lengthr;
				double nyl = (yl - yc) / lengthl; double nyr = (yr - yc) / lengthr;

				double dx = 0.1 * (fl * nxl + fr * nxr);
				double dy = 0.1 * (fl * nyl + fr * nyr);
				double dlength = std::hypot(dx, dy);
				if (err < dlength)  err = dlength;

				dxs[i] = dx; dys[i] = dy;

				i++;
			}

			std::size_t ii = 0;
			for (auto it = itbegin;	it != itend; it++)
			{				
				auto& [x, y, r] { ps[*it] };
				x += dxs[ii];		y += dys[ii];

				//output();
				double temp1 = (xb - xa) * (xb - x) + (yb - ya) * (yb - y);
				double temp2 = (xb - xa) * (x - xa) + (yb - ya) * (y - ya);
				double temp = std::pow(xb - xa, 2) + std::pow(yb - ya, 2);

				x = (temp1 * xa + temp2 * xb) / temp;
				y = (temp1 * ya + temp2 * yb) / temp;
				//output();
				ii++;
			}

			if (n % 50 == 0)
			{
				//output();				
				if(err  < ERR)
				break;
			}
		}	
	}
} 



void InitialMesh::meshGeneration(double extent, double ratio_Distance, double ratio_Radius, double ratio_Gradient)
{
	//Step 1: Initialization of point bubbles 

	updateConnectionsOfPoints();
	//output();
	
	setRadiiOfPointBubblesAccordingToRefinementExtent(extent);
	//output();

	updateRadiiOfPointBubblesAccordingToDistanceRatio(ratio_Distance);
	//output();

	updateRadiiOfPointBubblesAccordingToGradientRatio(ratio_Gradient);
	//output();

	
 //Step 2: Refinement of Boundary Edges 
 	double rmin = getPredefinedMinimumRadiusAccordingToRadiusRatio(ratio_Radius);
 
	for (std::size_t i = 0; i < 5; i++)
	{
		refineBoundaryEdgesAccordingToGradientRatio(ratio_Gradient);
		//output();

		updateConnectionsOfPoints();
		//output();

		updateRadiiOfPointBubblesAccordingToDistanceRatioAndMinimumRadius(ratio_Distance, rmin);
		//output();

		updateRadiiOfPointBubblesAccordingToGradientRatio(ratio_Gradient);
		//output();
	}
	
	//Step 3: Position Optimization of Local Bubbles 

	adjuestBoundaryPointsPositions();
	//output();

	//updateRadiiOfPointBubblesAccordingToGradientRatio(ratio_Gradient);
	//output();

	updateConnectionsOfPoints();
	//output();
}

