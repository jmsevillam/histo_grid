#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

typedef std::vector<double> vec;

void read_file(std::string filename,vec & x, vec & y,vec & z);
void init_vec(double & d,int N);
void dis(const vec & x,const  vec & y,const vec & z,const vec & x2,const  vec & y2,const vec & z2,vec & d);
double distance(double x,double y,double z,double x1,double y1,double z1);
double argmin(const vec & d);
void zero(vec & d);


void init_vec(vec & d,int N);

int main(int argc, char const *argv[]) {
    vec x,y,z,d;
    std::cout<<"creating vectors"<<std::endl;
    
    init_vec(x,atoi(argv[1]));
    init_vec(y,atoi(argv[1]));
    init_vec(z,atoi(argv[1]));

    std::cout<<"init vectors"<<std::endl;

    read_file("/all.dat",x,y,z);
    vec x2,y2,z2;
    init_vec(x2,atoi(argv[2]));
    init_vec(y2,atoi(argv[2]));
    init_vec(z2,atoi(argv[2]));
    init_vec(d,atoi(argv[2]));

    std::cout<<"init grid"<<std::endl;
    read_file("centers_6.dat",x2,y2,z2);
    std::cout<<"tessadsat "<<x[0]<<std::endl;

    dis(x,y,z,x2,y2,z2,d);

    return 0;
}

void read_file(std::string filename,vec & x, vec & y,vec & z){

    double xfile,yfile,zfile;

    std::stringstream ss;
    ss<<filename;
    std::fstream myfile(ss.str(), std::ios_base::in);

    for (int i = 0; !myfile.eof() ; i++){

        myfile >> xfile;
        myfile >> yfile;
        myfile >> zfile;

        x[i]=xfile;
        y[i]=yfile;
        z[i]=zfile;    

        if(i>=x.size()){
            myfile.close();
            break;
        }
    }
}

void init_vec(vec & dd,int N){
    dd.resize(N);
    zero(dd);
}

void zero(vec & dd){
for (int i = 0; i<dd.size();i++){
        dd[i]=0.0;        
    }
}

void dis(const vec & x,const  vec & y,const vec & z,const vec & x2,const  vec & y2,const vec & z2,vec & d){
    vec dd;
    init_vec(dd,x2.size());
        
    for (int i = 0; i<x.size();i++){
        if (i%1000==0){std::cout<<i<<std::endl;}
        zero(dd);
        for (int j = 0; j<x2.size();j++){
            dd[j]=distance(x[i],y[i],z[i],x2[j],y2[j],z2[j]);        
        }   
    d[argmin(dd)]+=1;  
    }

}

double distance(double x,double y,double z,double x1,double y1,double z1){
return sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1));
}


double argmin(const vec & d){
    double d_min=1.e10;
    int d_min_arg=-1;
    for (int i = 0; i<d.size();i++){
        if(d_min>d[i]){
            d_min=d[i];
            d_min_arg=i;
        }
    }
    return d_min_arg;
}
