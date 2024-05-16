#include <iostream>
#include <vector>
#include <stdexcept>
#include <tuple>

// Class representing points on the elliptic curve
class Point {
public:
    long long x, y;
    bool inf;

    Point(long long x = 0, long long y = 0, bool inf = false) : x(x), y(y), inf(inf) {}

    // Point at infinity
    static Point infinity() {
        return Point(0, 0, true);
    }

    bool operator==(const Point &other) const {
        return (inf && other.inf) || (x == other.x && y == other.y && inf == other.inf);
    }

    bool operator!=(const Point &other) const {
        return !(*this == other);
    }
};

// Class representing the elliptic curve y^2 = x^3 + ax + b over a finite field F_p
class EllipticCurve {
public:
    long long a, b, p;

    EllipticCurve(long long a, long long b, long long p) : a(a), b(b), p(p) {
        if (4 * a * a * a + 27 * b * b % p == 0) {
            throw std::invalid_argument("Invalid curve parameters");
        }
    }

    // Modular inverse using extended Euclidean algorithm
    long long mod_inv(long long n) const {
        long long t = 0, newt = 1;
        long long r = p, newr = n;
        while (newr != 0) {
            long long quotient = r / newr;
            std::tie(t, newt) = std::make_tuple(newt, t - quotient * newt);
            std::tie(r, newr) = std::make_tuple(newr, r - quotient * newr);
        }
        if (r > 1) throw std::runtime_error("No inverse exists");
        if (t < 0) t += p;
        return t;
    }

    // Point addition
    Point add(const Point &P, const Point &Q) const {
        if (P.inf) return Q;
        if (Q.inf) return P;
        if (P.x == Q.x && P.y != Q.y) return Point::infinity();
        
        long long m;
        if (P == Q) {
            m = (3 * P.x * P.x + a) * mod_inv(2 * P.y) % p;
        } else {
            m = (Q.y - P.y) * mod_inv(Q.x - P.x) % p;
        }
        
        long long x3 = (m * m - P.x - Q.x) % p;
        long long y3 = (m * (P.x - x3) - P.y) % p;
        return Point((x3 + p) % p, (y3 + p) % p);
    }

    // Scalar multiplication
    Point multiply(const Point &P, long long n) const {
        Point R = Point::infinity();
        Point Q = P;
        while (n > 0) {
            if (n % 2 == 1) {
                R = add(R, Q);
            }
            Q = add(Q, Q);
            n /= 2;
        }
        return R;
    }
};

int main() {
    // Example curve y^2 = x^3 + 2x + 3 over F_97
    EllipticCurve ec(2, 3, 97);
    Point P(3, 6);

    // Scalar multiplication
    long long k = 2;
    Point R = ec.multiply(P, k);
    std::cout << "2 * (" << P.x << ", " << P.y << ") = (" << R.x << ", " << R.y << ")\n";

    return 0;
}
