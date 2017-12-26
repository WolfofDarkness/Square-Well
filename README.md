# Square-Well
Calculates first order energy shifts using perturbation theory with two particles of same mass
//  file:shooting_method.cpp 
//
//  This program solves the schroedinger equation for an infinite square well using the shooting method(initially without normalization)
//
//  Programmer:  Andrew George george.494@osu.edu
//
//  Revision history:
//      29-05-2017: Original version 
//
//  Notes:  
//   * compile with: the corresponding make file, make_shooting_method
//

//*********************************************************************// 

// include files
#include <iostream> // this has the cout, cin definitions
#include <cmath> // useful math functions
#include <fstream> // allows me to stream things into data files
#include <iomanip> // use things such as setprecision and setw
using namespace std; // if omitted, then need std::cout, std::cin 

//*********************************************************************//


int
main ()
{
//first define the variables needed for finding the wave function and normalizing it. 
  double L = 1.0;
  double hc = 2.0;
  double psi = 0.0;
  double dpsi = 1.0;
  double ddpsi = 0.0;
  double E = 2.40;
  double x = 0;
  double dx = .01;
  double V = 0;






  ofstream out;
  out.open ("shooting_method.dat");
//now calculate the wave function 
  out << "E         x       psi" << endl;
//for (double E=0; E<=10; E +=.1)        //attempted to find a good E for it to work 
//{
  for (double i = 0.0; i <= L; i += .01)
    {
      x = 0; //for reinitialization 
      ddpsi = -2 * hc * (E - V) * psi; //second derivative of the wave function
      dpsi = dpsi + ddpsi * dx; //first derivative of the wave function 
      psi = psi + dpsi * dx; //wave function 
      x = x + i; //calculates the x value fo the particle from 0 to L

      out << fixed << setprecision (8) << setw (5) << E << " "
<< fixed << setprecision (8) << setw (5) << x << " "
<< fixed << setprecision (8) << setw (5) << psi << " " << endl;
    }
//}
  out.close ();

  return 0;
}

//  file:shooting_method_normalization.cpp 
//
//  This program solves the schroedinger equation for an infinite square well using the shooting method
//
//  Programmer:  Andrew George george.494@osu.edu
//
//
//  Notes:  
//   * compile with: the corresponding make file, make_shooting_method_normalization
//

//*********************************************************************// 

// include files
#include <iostream> // this has the cout, cin definitions
#include <cmath> // useful math functions
#include <fstream> // allows me to stream things into data files
#include <iomanip> // use things such as setprecision and setw
using namespace std; // if omitted, then need std::cout, std::cin 

//*********************************************************************//


int
main ()
{
//first define the variables needed for finding the wave function and normalizing it. 
  double L = 1.0;                  //length of well
  double hc = 2.0;                 
  double psi = 0.0;                //guess wavefunction
  double dpsi = 1.0;               //first derivative of the wavefunction
  double ddpsi = 0.0;              //second derviative of the wavefunction
  double E = 2.40;                 //Energy found through shooting method
  double x = 0;
  double dx = .01;
  double V = 0;                     //potential energy
  double xt = 0;
  double b = 4.37;                  //test normalization factor for numerical method. Change around to find correct factor. 
  double wavefunct = 0;             //known wavefunction
  double normalize = 0;
  double A = 0;
  double pi=atan(1.0)*4.0;
  double psinormal=0;

//will first integrate the analytical solution to find a normalized probability density
  ofstream out2;
out2.open ("analytical_solution.dat"); 
out2 << "x*t       wavefunction^2"<<endl;
for (double j=0.0; j<=L; j +=.01)
{
xt=0; 
wavefunct=sqrt(2/L)*sin(pi*j/L);
normalize=normalize+(wavefunct*wavefunct*dx);
xt=xt+j;

out2<<fixed<<setprecision(8)<<setw(5)<<xt<< " "
    <<fixed<<setprecision(8)<<setw(5)<<wavefunct*wavefunct<< " "
     <<endl; 
}
out2.close(); 

cout<<" The integrated probability density of the analytical function is:"<<normalize<<endl;




//now will integrate my numerical algorithm to find to normalize it. 
ofstream out; 
out.open ("numerical_solution.dat"); 
out<< " x          psinormal^2 "<<endl; 
for (double i=0.0; i<=L; i +=.01) 
{
 x=0;                      //for reinitialization 
 ddpsi=-2*hc*(E-V)*psi;    //second derivative of the wave function
 dpsi=dpsi+ddpsi*dx;        //first derivative of the wave function 
 psi=psi+dpsi*dx;           //wave function 
 x=x+i;                    //calculates the x value fo the particle from 0 to L

 psinormal=b*psi; 
 
 A=A+psinormal*psinormal*dx;
 

out<<fixed<<setprecision(8)<<setw(5)<<x<< " "
   <<fixed<<setprecision(8)<<setw(5)<<psinormal*psinormal<< " "
   <<endl; 
} 
out.close(); 

cout<<" The integrated probability density of the numerical function is:"<<A<<endl;
return 0; 
}

//  file:first_order_petrub.cpp
//
//  This program solves the schroedinger equation for an infinite square well using the shooting method
//
//  Programmer:  Andrew George george.494@osu.edu
//
//
//  Notes:  
//   * compile with: the corresponding make file, make_first_order_petrub
//

//*********************************************************************// 

// include files
#include <iostream> // this has the cout, cin definitions
#include <cmath> // useful math functions
#include <fstream> // allows me to stream things into data files
#include <iomanip> // use things such as setprecision and setw
using namespace std; // if omitted, then need std::cout, std::cin 

//*********************************************************************//


int
main ()
{
int scalar;

cout<<"For this first order energy shift, we will be using the values L=1.0, hc=2.0, psi=0.0, dpsi=1.0, ddpsi=0.0, x=0, dx=.01, V=0, E=2.40, and b=4.37"<<endl;
cout<< "The energy and normalization factor were found through the shooting methods. Please look at shooting_method.cpp and shooting_method_normalization.cpp for more details."<<endl;
cout<<endl;
cout<<endl;

cout<<"This will program calculate first order energy shifts for a single particle with a pertubative Hamiltonian that is a delta function and has a scalar factor. Please input the scalar factor."<<endl;


cin>>scalar; 

cout<<"Your first order energy shifts are:"<<(2*scalar/1.0)<<endl;
cout<<"Your total energy is:"<<(2*scalar/1.0)+2.40<<endl; 

//Now for symmetrization


return 0; 

}
