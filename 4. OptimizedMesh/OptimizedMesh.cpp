
#include "OptimizedMesh.h"
//#include <numbers>
//using std::numbers::pi;

// mesh data input 
void OptimizedMesh::input()
{
	std::ifstream ist;
	ist.open("Mesh.txt");
	std::string ss;	

	// input Points
	std::getline(ist, ss);
	size_t np{ 0 }; ist >> np;
	for (size_t i = 0; i < np; i++)
	{
		static size_t n{ 0 }; 	 
		static double x{ 0.0 }, y{ 0.0 }, r{0.0};
		ist >> n >> x >> y >> r;
		ps.emplace_back(x, y, r);
	}
	ist >> std::ws;


	// input Triangles 
	std::getline(ist, ss);
	size_t nt{ 0 }; ist >> nt;	
	for (size_t i = 0; i < nt; i++)
	{
		static size_t n{ 0 };
		static size_t a{ 0 }, b{ 0 }, c{ 0 };
		ist >> n >> a >> b >> c;
		ts.emplace_back(i, a, b, c);
	}
	ist >> std::ws;

	ist.close();
}


void OptimizedMesh::output() const
{
	const std::string name{ "OptimizedMesh.txt" };
	std::ofstream ost;
	ost.open(name);
	// output Points 
	ost << "Points\n";
	ost << ps.size() << '\n';
	for (size_t i = 1; auto & p : ps)
	{
		ost << std::scientific << i << '\t'
			<< p.x << '\t' << p.y << '\t' << p.r << '\n';
		i++;
	}

	// output Triangles 
	ost << "Triangles\n";
	ost << ts.size() << '\n';
	for (size_t i = 1; auto & t : ts)
	{
		ost << std::scientific << i << '\t'
			<< t.a << '\t' << t.b << '\t' << t.c << '\t' << '\n';
		i++;
	}

	ost.close();
}

void OptimizedMesh::eraseTriangeEdges(std::size_t n)
{
	for (auto& [e, ep] : ts[n].eep3)
	{
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


void OptimizedMesh::insertTriangeEdges(std::size_t n)
{
	for (auto& [e, ep] : ts[n].eep3)
		eepss[e].emplace_back(ep);
}


//data_initialization
void OptimizedMesh::initialization()
{
	//eepss initialization
	for (std::size_t i=0; auto& t : ts)
	{
		insertTriangeEdges(i);
		i++;
	}

		
	//Points initialization (set radii) 
	for (auto& p : ps)	
		p.r *= 1.7;	

	for (auto& [e, eps] : eepss)
	{
		if (eps.size() == 1)
		{
			ps[e.a].isatboundary = true;
			ps[e.b].isatboundary = true;
		}
	}




	//Neighbouring Points
	std::vector<std::unordered_set<size_t>> apss(ps.size());
		
	//Points initialization (set nbs and aps )
	for (const auto& t : ts)
	{
		for (auto& [e, ep] : t.eep3)
		{
			auto& [a, b] {e};

			apss[a].insert(b);		
			apss[b].insert(a);	
		}
	}


	//Points initialization (set aaps) 
	for (std::size_t i = 0; auto & p:ps)
	{	
		p.aaps = apss[i];
		i++;
	}
		

	std::size_t LEVEL{ 3 }; // Give 3 or 4
	for (std::size_t i = 1; i < LEVEL; i++)
	{
		for (std::size_t ii = 0; auto& p : ps)
		{			
			auto& aaps{ p.aaps }; auto aaps_temp{ aaps };
			for (auto aa:aaps_temp)
			{		
				for (auto& a : apss[aa])					
					aaps.insert(a);
			}	
			ii++;
		}
	}

	for (size_t i = 0; auto& p : ps)
	{
		auto& aaps{ p.aaps };	
		aaps.erase(i);		
		i++;
	}

	output();
}


//new positions calculation
void OptimizedMesh::new_positions_calculation()
{
	for (auto& p : ps)
	{
		if (p.isatboundary == false)
		{
			auto& dx = p.dx; dx = 0.0;
			auto& dy = p.dy; dy = 0.0;
			for (auto& a : p.aaps)
			{
				auto& pa{ ps[a] };
				double length = std::hypot(pa.x - p.x, pa.y - p.y);
				double r = p.r + pa.r;
				
				double f = (length > r) ? 0.0 : (length - r);
				double nx = (pa.x - p.x) / length;
				double ny = (pa.y - p.y) / length;
				dx += 0.1 * f * nx;
				dy += 0.1 * f * ny;
			}			
		}
	}

	for (auto& p : ps)
	{
		if (p.isatboundary == false)
		{
			p.x += p.dx;
			p.y += p.dy;
		}		
	}
}


//radii calculation
void OptimizedMesh::radii_calculation()
{
	for (auto&& p : ps)
	{
		if (p.isatboundary == false)
		{		
			double sum_r{ 0.0 }, sum{ 0.0 };
			for (auto& a : p.aaps)
			{
				auto& pa{ ps[a] };
				auto length = std::hypot(pa.x - p.x, pa.y - p.y);
				double r = p.r + pa.r;
				if (length < r)
				{
					sum_r += pa.r / length;
					sum += 1.0 / length;
				}
			}

			p.dr = sum_r / sum - p.r;
		}
	}

	for (auto&& p : ps)
	{
		if (p.isatboundary==false)			
		p.r+= p.dr;
	}
}

double OptimizedMesh::getRelativeLength(std::size_t a, std::size_t b, std::size_t c1, std::size_t c2) const
{
	auto& pa{ ps[a] }; auto& pb{ ps[b] };
	auto& pc1{ ps[c1] }; auto& pc2{ ps[c2] };

	double xa{ pa.x }, xb{ pb.x }, ra{pa.r}, ya{ pa.y }, yb{ pb.y }, rb{ pb.r };
	double xc1{ pc1.x }, yc1{ pc1.y }, rc1{ pc1.r }, xc2{ pc2.x }, yc2{ pc2.y }, rc2{ pc2.r };

	double rlab =  std::hypot(xa - xb, ya - yb) / (ra + rb);

	double rlc1c2 =  std::hypot(xc1 - xc2, yc1 - yc2) / (rc1 + rc2)  ;
		
	return rlc1c2 - rlab;

}



void OptimizedMesh::updatePointsConnections()
{
	using RelativeLengthEs = std::multimap<double, const Edge>;
	using EIts = std::unordered_map<Edge, RelativeLengthEs::iterator, EdgeHash, EdgeEqual>;
	RelativeLengthEs rles;  EIts eits;

	for (auto& [e, eps] : eepss)
	{
		if (eps.size() == 2)
		{
			double rl= getRelativeLength(e.a, e.b, eps.front().c, eps.back().c);
			if (rl < - 1e-8)
			{
				auto it = rles.emplace(rl, e);
				eits.emplace(e, it);
			}
		}
	}

	size_t i = 0;
	while (rles.size() > 0)
	{
		auto itmin = rles.begin();
		auto& emin{ itmin->second };	auto& epsmin = eepss.at(emin);
		auto [a, b] {emin};
		auto [nt1, c1] { epsmin.front() };
		auto [nt2, c2] { epsmin.back()};


		//delete olds 
		for (auto n : { nt1,nt2 })
			eraseTriangeEdges(n);


		Edges eolds;
		for (auto n : { nt1,nt2 })
			for (auto& [e, ep] : ts[n].eep3)
				eolds.emplace(e);

		for (auto& e : eolds)
		{
			auto it = eits.find(e);
			if (it != eits.end())
			{
				rles.erase(it->second);
				eits.erase(it);
			}
		}


		//add news
		ts[nt1] = Triangle(nt1, c1, c2, a);
		ts[nt2] = Triangle(nt2, c1, c2, b);

		for (auto n : { nt1,nt2 })
			insertTriangeEdges(n);


		Edges enews;
		for (auto n : { nt1,nt2 })
			for (auto& [e, ep] : ts[n].eep3)
				enews.emplace(e);

		for (auto& e : enews)
		{
			auto& eps = eepss.at(e);
			if (eps.size() == 2)
			{
				double rl = getRelativeLength(e.a, e.b, eps.front().c, eps.back().c);
				if (rl <  - 1e-8)
				{
					auto it = rles.emplace(rl, e);
					eits.emplace(e, it);
				}
			}
		}

		//if ((i++) % (20) == 0)
		//{
		//	output();
		//	continue;
		//}

	}
}





void OptimizedMesh::optimization()
{
	//double r_max{ 0.0 };

	//size_t NN {10 * static_cast<size_t> (std::sqrt(ps.size())) };
	size_t NN{10* static_cast<size_t> (std::sqrt(ps.size())) };
	for (size_t i = 0; i < NN; i++)
	{
		for (size_t j = 0; j < 20; j++)
			new_positions_calculation();

		radii_calculation();		
		
		//if (i % 10 == 0)
		//{
		//	output();
		//	continue;
		//}

	}

	//output();
	updatePointsConnections();

}

