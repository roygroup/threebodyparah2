#ifndef _JPOINT_H
#define _JPOINT_H

#include <tuple>

class JacobiPoint {
private:
	double R, s, u;
	
	JacobiPoint(double R, double s, double u);
	
public:
	// constructor and destructor
	static JacobiPoint from_jacobicoords(double, double, double);
	static JacobiPoint from_unsorted_pairdistances(double, double, double);
	
	JacobiPoint(const JacobiPoint&);
	JacobiPoint& operator=(const JacobiPoint&);
	~JacobiPoint();
	
	double get_R() const;
	double get_s() const;
	double get_u() const;

    std::tuple<double, double, double> unpack() const;
};

#endif
