//
//  main.cpp
//  Assignment_#1
//
//  Created by Isra Ragheb on 10/3/18.
//  Copyright © 2018 DBB. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <cmath>

struct point{
    double x;
    double y;
    double z;
};

//Assignment1 (Part a)
 std::tuple<struct point,struct point> getSphereintersect(struct point center,double radius,struct point pnode, struct point pfinal,int method){
   struct point PIntersection1;
   struct point PIntersection2;
   struct point v;
    v.x=pfinal.x-pnode.x;
    v.y=pfinal.y-pnode.y;
    v.z=pfinal.z-pnode.z;
    if(method==2){
        //Algebric method
        
        //||v||2
        double a=sqrt(v.x*v.x+v.y*v.y+v.z*v.z)*sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
        struct point vb;
        struct point bray;
        
        //2v.[P0-O]
        vb.x=2*v.x;
        vb.y=2*v.y;
        vb.z=2*v.z;
        bray.x=pnode.x-center.x;
        bray.y=pnode.y-center.y;
        bray.z=pnode.z-center.z;
        double b=(vb.x*bray.x)+(vb.y*bray.y)+(vb.z*bray.z);
        
        //[P0-O]2-r2
        double c=((bray.x*bray.x)+(bray.y+bray.y)+(bray.z*bray.z))-(radius*radius);
        
        //Get the value of t
        double delta=(b*b)-4*a*c;
        if(delta<0){
            PIntersection1.x=-1;
            PIntersection1.y=-1;
            PIntersection1.z=-1;
        }
        else{
            if(delta==0){
                 double tPlus=((-1*b))/(2*a);
                PIntersection1.x=pnode.x+(tPlus*v.x);
                PIntersection1.y=pnode.y+(tPlus*v.y);
                PIntersection1.z=pnode.z+(tPlus*v.z);
                
                PIntersection2.x=0;
                PIntersection2.y=0;
                PIntersection2.z=0;
            }
            else{
                double tminus=((-1*b)-(sqrt((b*b)-(4*a*c))))/(2*a);
                double tPlus=((-1*b)+(sqrt((b*b)-(4*a*c))))/(2*a);
                
                PIntersection1.x=pnode.x+(tPlus*v.x);
                PIntersection1.y=pnode.y+(tPlus*v.y);
                PIntersection1.z=pnode.z+(tPlus*v.z);
                
                PIntersection2.x=pnode.x+(tminus*v.x);
                PIntersection2.y=pnode.y+(tminus*v.y);
                PIntersection2.z=pnode.z+(tminus*v.z);

            }
        }
        
        
    }else{
        if(method==1){
            //Geomatric
            //a=center-pnode
            struct point a;
            a.x=center.x-pnode.x;
            a.y=center.y-pnode.y;
            a.z=center.z-pnode.z;
            
            //standerizing v
            double magv=sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
            struct point standv;
            standv.x=(1/magv)*v.x;
            standv.y=(1/magv)*v.y;
            standv.z=(1/magv)*v.z;
            
            //appling geo formula
            double aa=(a.x*a.x)+(a.y*a.y)+(a.z*a.z);
            double av=(a.x*standv.x)+(a.y*standv.y)+(a.z*standv.z);
            double cons=(av)-sqrt((radius*radius)-((aa)-(av*av)));
            PIntersection1.x=pnode.x+(cons*standv.x);
            PIntersection1.y=pnode.y+(cons*standv.y);
            PIntersection1.z=pnode.z+(cons*standv.z);
            
            PIntersection2.x=0;
            PIntersection2.y=0;
            PIntersection2.z=0;
        }
        else{
            PIntersection1.x=0;
            PIntersection1.y=0;
            PIntersection1.z=0;
            
            PIntersection2.x=0;
            PIntersection2.y=0;
            PIntersection2.z=0;
            
        }
    }
    
    return std:: make_tuple(PIntersection1,PIntersection2);
}


//Assignment1 (Part b)
std::tuple<struct point,struct point> getConeIntersection(struct point centerC,double radius,double heightC,struct point pnode,struct point p2){
    
    struct point PInescetion;
    struct point PInescetion1;
    
    //normalizatoin the D vector
    struct point unitD, D;
    unitD.x=p2.x-pnode.x;
    unitD.y=p2.y-pnode.y;
    unitD.z=p2.z-pnode.z;
    
    double mag=sqrt(pow(unitD.x, 2)+pow(unitD.y, 2)+pow(unitD.z, 2));
    D.x=(1/mag)*unitD.x;
    D.y=(1/mag)*unitD.y;
    D.z=(1/mag)*unitD.z;
    
    //Geting the point C(the tip of the cone)
    struct point C;
    C.x=centerC.x;
    C.y=centerC.y+ heightC;
    C.z=centerC.z;
    
    //normalizatoin the V vector (Center-Cpont)
    struct point unitv,v ;
    unitv.x=centerC.x-C.x;
    unitv.y=centerC.y-C.y;
    unitv.z=centerC.z-C.z;
    double magv=sqrt(pow(unitv.x, 2)+pow(unitv.y, 2)+pow(unitv.z, 2));
    v.x=(1/magv)*unitv.x;
    v.y=(1/magv)*unitv.y;
    v.z=(1/magv)*unitv.z;
    
    
    
    struct point co;
    co.x=pnode.x-C.x;
    co.y=pnode.y-C.y;
    co.z=pnode.z-C.z;
    
    double theta=atan(radius/heightC)*(180/3.142);
    
    double dv=(D.x*v.x)+(D.y*v.y)+(D.z*v.z);
    //calculating A
    double a=(dv*dv)-(cos(theta)*cos(theta));
    double cov=(co.x*v.x)+(co.y*v.y)+(co.z*v.z);
    double cod=(co.x*D.x)+(co.y*D.y)+(co.z*D.z);
    double cdv=(C.x*dv)+(C.y*dv)+(C.z*dv);
    //calculating B
    double b=2*((cdv*cov)-cod*cos(theta)*cos(theta));
    //calculating C
    double c=(cov*cov)-((co.x*co.x)+(co.y*co.y)+(co.z*co.z))*(cos(theta)*cos(theta));
    
    
    double delta=(b*b)-4*a*c;
    
    //Theta should be less that 90 because we squared the equation
    if(theta<90){
        //checking on delta to find out the number of possiable intersections
        if(delta<0){
            //no intersectin
            PInescetion.x=0;
            PInescetion.y=0;
            PInescetion.z=0;
        }
        else{
            if(delta==0){
                //single intersection
                double t=(-b*a)/(2*a);
                PInescetion.x=pnode.x+t*D.x;
                PInescetion.y=pnode.y+t*D.y;
                PInescetion.z=pnode.z+t*D.z;
            }
            else{
                //to possiable intersction
                double tminus=((-1*b)-(sqrt((b*b)-(4*a*c))))/(2*a);
                double tPlus=((-1*b)+(sqrt((b*b)-(4*a*c))))/(2*a);
                PInescetion.x=pnode.x+tminus*D.x;
                PInescetion.y=pnode.y+tminus*D.y;
                PInescetion.z=pnode.z+tminus*D.z;
                
                PInescetion1.x=pnode.x+tPlus*D.x;
                PInescetion1.y=pnode.y+tPlus*D.y;
                PInescetion1.z=pnode.z+tPlus*D.z;
                
            }
        }
    }
    
   
    return std::make_tuple(PInescetion,PInescetion1);
}



int main(int argc, const char * argv[]) {
    char method;
    struct point intr1, intr2;
    double radius=0;
    double height=0;
    int spheremethod=1;
    std::cout <<"Please enter which method you wish to test(Method A or B)";
    std::cin >> method;
    struct point center;
    struct point node;
    struct point end;
    
    if (method=='A') {
        //Getting the center of the sphere from the user
        std::cout << "Please enter the center of the sphere:" << std::endl;
        std::cin >> center.x;
        std::cin >> center.y;
        std::cin >> center.z;
        
        //Getting the start of the tracing ray from the user
        std::cout << "Please enter the start of the tracing ray:" << std::endl;
        std::cin >> node.x;
        std::cin >> node.y;
        std::cin >> node.z;
        
        //Getting the end of the tracing ray from the user
        std::cout << "Please enter the end of the tracing ray:" << std::endl;
        std::cin >> end.x;
        std::cin >> end.y;
        std::cin >> end.z;
        
        //Getting the center of the ray from the user
        std::cout << "Please enter the radius of the sphere:" << std::endl;
        std::cin >> radius;
        
        //Getting the method that they wish to use from the user
        std::cout << "Please enter the method you wish to use(1 for the geometric method and 2 for for the algebric method):" << std::endl;
        std::cin >> spheremethod;
        
        std::cout <<std:: endl;
        
        //calling the method for calculating the ray sphere intersetion
        std:: tie(intr1,intr2)=getSphereintersect(center, radius, node, end, spheremethod);
        
    }
    else{
        if(method=='B'){
            //Getting the center of the base of cone the from the user
            std::cout << "Please enter the center of the base of the cone:" << std::endl;
            std::cin >> center.x;
            std::cin >> center.y;
            std::cin >> center.z;
            
            //Getting the start of the tracing ray from the user
            std::cout << "Please enter the start of the tracing ray:" << std::endl;
            std::cin >> node.x;
            std::cin >> node.y;
            std::cin >> node.z;
            
            //Getting the end of the tracing ray from the user
            std::cout << "Please enter the end of the tracing ray:" << std::endl;
            std::cin >> end.x;
            std::cin >> end.y;
            std::cin >> end.z;
            
            //Getting the raduis of the cone from the user
            std::cout << "Please enter the radius of the cone:" << std::endl;
            std::cin >> radius;
            
            //Getting the height of the cone from the user
            std::cout << "Please enter the height of the cone:" << std::endl;
            std::cin >> height;
            
             std::cout <<std:: endl;
            
            //calling the method for calculating the ray cone intersetion
            std:: tie(intr1,intr2)=getConeIntersection(center, radius, height, node, end);
        }
    }
    
    //Priting the intersection poin to the user
    std::cout << intr1.x<< std:: endl;
    std::cout << intr1.y<< std:: endl;
    std::cout << intr1.z<<std:: endl;
    
    std::cout <<std:: endl;
    
    std::cout << intr2.x<< std:: endl;
    std::cout << intr2.y<< std:: endl;
    std::cout << intr2.z<<std:: endl;
    
    
    std::cout << "Hello, World!\n";
    
    return 0;
}
