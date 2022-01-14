#include <cmath>
#include <tuple>

#include "utils.h"

namespace utils
{
    double cosine_transition(double x, double xmin, double xmax) {
        if      (x <= xmin) {
            return 0.0;
        }
        else if (x >= xmax) {
            return 1.0;
        }
        else {
            double k = (x - xmin)/(xmax - xmin);
            return 0.5 * (1.0 - cos(M_PI*k));
        }
    }
}

namespace atm
{
    double wterm(const JacobiPoint& jpoint) {
        double s = jpoint.get_s();
        double u = jpoint.get_u();

        return wterm(s, u);
    }

    double wterm(double s, double u) {
        double cosu   = cos(HALFPI*u);
        double w_sqrt = s * (cosu + sqrt(3.0 + cosu*cosu));

        return w_sqrt*w_sqrt;
    }

    double tterm(const JacobiPoint& jpoint) {
        double s = jpoint.get_s();
        double u = jpoint.get_u();

        return tterm(s, u);
    }

    double tterm(double s, double u) {
        double cos2u = cos(2.0*HALFPI*u);
        double w     = wterm(s, u);

        return 1.0 - 2.0*w*cos2u + w*w;
    }

    double fterm(const JacobiPoint& jpoint) {
        double u = jpoint.get_u();

        double w = wterm(jpoint);
        double t = tterm(jpoint);

        double cosu  = cos(HALFPI*u);
        double numerator = (1.0 - w*cosu*cosu) * (1.0 - w);

        double f = -3.0 * numerator/t;

        return f;
    }

    double norm(const JacobiPoint& jpoint) {
        double R = jpoint.get_R();
        double t = tterm(jpoint);
        
        return (pow(R, 9) * pow(t, 1.5)) / (64.0 * C9);
    }

    double invnorm(const JacobiPoint& jpoint) {
        double R = jpoint.get_R();
        double t = tterm(jpoint);
        
        return (64.0 * C9) / (pow(R, 9) * pow(t, 1.5));
    }

    double energy(const JacobiPoint& jpoint) {
        return invnorm(jpoint) * (1.0 + fterm(jpoint));
    }
}
