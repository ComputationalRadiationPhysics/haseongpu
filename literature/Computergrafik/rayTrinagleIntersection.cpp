/* 
 * File:   main.cpp
 * Author: Frank Liebold
 *
 * Created on 6. Juni 2012, 11:13
 */




#include <cstdlib>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

/*
 * 
 */

//aufbauend auf
//Fast, Minimum Storage Ray/Triangle Intersection
//von Tomas Moeller und Ben Trumbore
//siehe Google


//Vektorstruktur
struct Vek3d {
private:
    //Membervariablen
    double _x;
    double _y;
    double _z;

public:

    //Konstruktor
    //inline Vek3d(double X=0.0f,double Y=0.0f, double Z=0.0f):_x(X),_y(Y),_z(Z){}

    Vek3d(double const &X, double const &Y, double const &Z) : _x(X), _y(Y), _z(Z) {
    }

    Vek3d() : _x(0.), _y(0.), _z(0.) {
    }

    //Destruktor

    ~Vek3d() {
    }

    //Kopierkonstruktor

    Vek3d(const Vek3d &v) {
        _x = v.x();
        _y = v.y();
        _z = v.z();
    }

    ////////////////////////////////////////////
    //Operatoren

    //Zuweisungsoperator

    const Vek3d& operator=(Vek3d const &v) {
        _x = v.x();
        _y = v.y();
        _z = v.z();
        return *this;
    }

    //Vergleichsoperator

    inline bool operator==(Vek3d const &v) {
        if (_x == v.x() && _y == v.y() && _z == v.z())
            return 1;
        else return 0;
    }

    /** Subtraktionsoperator */
    inline Vek3d operator-(Vek3d const &v) const {
        return Vek3d(_x - v.x(), _y - v.y(), _z - v.z());
    }

    inline Vek3d& operator-=(Vek3d const &v) {
        _x -= v.x();
        _y -= v.y();
        _z -= v.z();
        return *this;
    }

    // Additionsoperator

    inline Vek3d operator+(Vek3d const &v) const {
        return Vek3d(_x + v.x(), _y + v.y(), _z + v.z());
    }

    inline Vek3d& operator+=(Vek3d const &v) {
        _x += v.x();
        _y += v.y();
        _z += v.z();
        return *this;
    }

    // Multiplikationsoperator mit Skalar

    inline Vek3d operator*(double const &s) const {
        return Vek3d(s*_x, s*_y, s * _z);
    }

    //Multiplikation mit Skalar

    inline Vek3d& operator*=(double const &d) {
        _x *= d;
        _y *= d;
        _z *= d;
        return *this;
    }


    // Divisionsoperator Division durch Skalar

    inline Vek3d operator/(double const &s) const {
        return operator*(1.0f / s);
    }

    inline Vek3d& operator/=(double const &d) {
        return *this *= 1.0f / d;
    }

    //Skalarprodukt

    inline double operator*(Vek3d const &v) const {
        return v.x() * _x + v.y() * _y + v.z() * _z;
    }

    //Kreuzprodukt

    inline Vek3d operator^(Vek3d const &v) const {
        return Vek3d(_y * v.z() - _z * v.y(), _z * v.x() - _x * v.z(), _x * v.y() - _y * v.x());
    }

    // Output Stream

    friend std::ostream& operator<<(std::ostream &os, const Vek3d& v) // Static?!
    {
        return os << '[' << v.x() << ", " << v.y() << ", " << v.z() << ']';
    }

    inline const double& x() const {
        return _x;
    }

    inline double& x() {
        return _x;
    }

    inline const double& y() const {
        return _y;
    }

    inline double& y() {
        return _y;
    }

    inline const double& z() const {
        return _z;
    }

    inline double& z() {
        return _z;
    }
};




//###################################################################
//Schnitt Gerade Dreieck
//Input:
//Schnittpunkt-Vektor schnittp -> wird berechnet
//Geraden-Ortsvektor (Ursprung) rayOrigin
//Geraden-Richtungsvektor rayDirection
//Dreieckspunktvektoren: p1,p2,p3
//barizentrische Koordinaten des Schnittpunktes s2,s3 (p=s1*p1+s2*p2+s3*p3)-> werden berechnet
//s1 waere s1=1-s2-s3
//t: Vielfaches des Richtungsvektors der Geraden
//Output:
//1 wenn Schnittpunkt vorhanden und innerhalb Dreieck,ansonsten 0

bool intersectionRayTriangle(Vek3d &intersectionPoint, //Schnittpunkt
        const Vek3d &rayOrigin, //Ursprung des Strahls
        const Vek3d &rayDirection, //Richtungsvektor des Strahls
        const Vek3d &p1, //1.Punkt des Dreiecks
        const Vek3d &p2, //2.Punkt des Dreiecks
        const Vek3d &p3, //3.Punkt des Dreiecks
        double& s2, //2.barizentrische Koordinate des Dreiecks
        double& s3, //3.barizentrische Koordinate des Dreiecks
        //1.barizentrische Koordinate des Dreiecks ergibt sich mit 1.-s2-s3
        double &t, //Geradenparameter
        bool berechneSchnittp = 1)//soll Schnittpunkt berechnet werden (1) oder nicht (0)
{
    double eps = 1e-8;
    double determinante;

    //side12 und side13 sind Vektoren der Seiten
    //cross ist eine Hilfsvariable
    Vek3d side12, side13, cross;

    //Berechnung von Vektoren der Seiten:
    //1.Seite side12 von p1 nach p2
    side12.x() = p2.x() - p1.x();
    side12.y() = p2.y() - p1.y();
    side12.z() = p2.z() - p1.z();

    //2.Seite side13 von p1 nach p3
    side13.x() = p3.x() - p1.x();
    side13.y() = p3.y() - p1.y();
    side13.z() = p3.z() - p1.z();

    //Gleichsetzen von Gereadengleichung und Ebenengleichung
    //Modell:  sp=p1+s2*(p2-p1)+s3*(p3-p1)=rayOrigin+t*rayDirection
    //Berechnung mit Cramerscher Regel zum Loesen des Gleichungssystems:
    //s2*(p2-p1)+s3*(p3-p1)-t*rayDirection = rayOrigin-p1 -> zu bestimmende Parameter: s2,s3,-t 

    //Kreuzprodukt von side13 und rayDirection
    cross.x() = side13.y() * rayDirection.z() - side13.z() * rayDirection.y();
    cross.y() = side13.z() * rayDirection.x() - side13.x() * rayDirection.z();
    cross.z() = side13.x() * rayDirection.y() - side13.y() * rayDirection.x();

    //Berechnung der Determinante mit Skalarprodukt
    determinante = cross.x() * side12.x() + cross.y() * side12.y() + cross.z() * side12.z();

    if (determinante>-eps && determinante < eps) return 0;

    //Abstand Ursprung des Strahls zu p1
    Vek3d p1_rayOrigin; //=rayOrigin-p1;
    p1_rayOrigin.x() = rayOrigin.x() - p1.x();
    p1_rayOrigin.y() = rayOrigin.y() - p1.y();
    p1_rayOrigin.z() = rayOrigin.z() - p1.z();

    //barizentrische Koordinaten
    // sp=s1*p1+s2*p2+s3*p3
    //2. barizentrische Koordinate s2
    //=Skalarprodukt von p1_rayOrigin und cross
    s2 = cross.x() * p1_rayOrigin.x() + cross.y() * p1_rayOrigin.y() + cross.z() * p1_rayOrigin.z();

    //Hilfsvariable
    Vek3d tempcross;
    //zunaenaehst Kreuzprodukt von rayDirection und side12
    tempcross.x() = rayDirection.y() * side12.z() - rayDirection.z() * side12.y();
    tempcross.y() = rayDirection.z() * side12.x() - rayDirection.x() * side12.z();
    tempcross.z() = rayDirection.x() * side12.y() - rayDirection.y() * side12.x();

    //s3=Skalarprodukt von rayDirection und side12
    //s3=(rayDirection x side12) *p1_rayOrigin
    s3 = tempcross.x() * p1_rayOrigin.x() + tempcross.y() * p1_rayOrigin.y() + tempcross.z() * p1_rayOrigin.z();

    //Cramersche Regel -> Division durchfuehren
    double invdet = 1. / determinante;

    s2 = invdet*s2;
    s3 = invdet*s3;

    //weitere Verwendung der Hilfsvariable fuer Berechnung des Geradenparameters t
    //zunaechst Kreuzprodukt von side13 und side12
    tempcross.x() = side13.y() * side12.z() - side13.z() * side12.y();
    tempcross.y() = side13.z() * side12.x() - side13.x() * side12.z();
    tempcross.z() = side13.x() * side12.y() - side13.y() * side12.x();

    //t ist dann das Skalarprodukt von tempcross und p1_rayOrigin
    //t=(seite13,seite12) *p1_rayOrigin = -(seite12,seite13) *p1_rayOrigin
    t = tempcross.x() * p1_rayOrigin.x() + tempcross.y() * p1_rayOrigin.y() + tempcross.z() * p1_rayOrigin.z();

    t = invdet*t;

    //Test,ob der Schnittpunkt innerhalb des Dreiecks liegt:

    //Ueberschereitungstest fuer barizentrische Koordinaten
    if (s2 < 0. || s2 > 1.) return 0;

    //Ueberschereitungstest fuer barizentrische Koordinaten
    if (s3 < 0. || s3 > 1.) return 0;

    //0 <= s1=1-s2-s3 <= 1 -> s2+s3<1   (s2+s3>0 schon durchgefuehrt,da s2>0 s3>0)
    if (s2 + s3 > 1.) return 0;

    //Test, ob Strahl in Richtung des Dreiecks zeigt:
    if (t < 0.) return 0;

    //Schnittpunktberechnung
    if (berechneSchnittp) {
        intersectionPoint.x() = rayOrigin.x() + t * rayDirection.x();
        intersectionPoint.y() = rayOrigin.y() + t * rayDirection.y();
        intersectionPoint.z() = rayOrigin.z() + t * rayDirection.z();
    }

    return 1;
}




//#######################################################
int main(int argc, char** argv) {

    //Strahl
    Vek3d rayDirection(1., 0.1, 0.); //rayDirection(1,0,0)
    Vek3d rayOrigin(-4., 0., -0.3); //rayOrigin(-4,0,0)
    
    //Dreieckspunkte
    Vek3d p1(1., -1., -1.);
    Vek3d p2(1., 1., -1.);
    Vek3d p3(1., 0., 1.);

    //Schnittpunkt
    Vek3d sp;
    //Parameter
    //s2,s3 baryzentrische Koordinaten des Dreiecks
    //t Geradenparameter (Strahlmodell)
    double s2, s3, t;
    
    
    if (intersectionRayTriangle(sp, rayOrigin, rayDirection, p1, p2, p3, s2, s3, t)) {
        cout << "sp=" << sp << "\ns1=" << 1.- s2 - s3 << "   s2=" << s2 
             << "   s3=" << s3 << "   t=" << t << endl;

        cout << "\nTest (sollte Schnittpunkt entsprechen):\nhier berechnet mit baryzentrischen Koordinaten\n" 
             << p1 * (1. - s2 - s3) + p2 * s2 + p3 * s3 << endl;

    } else cout << "fehlgeschlagen, kein Schnittpunkt\n";


    return 0;
}

