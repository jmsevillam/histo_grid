#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include <mpi.h>

typedef std::vector<double> vec;

void read_file(std::string filename,double* x, double* y,double* z,int, int, int);
void dis(const double* x,const double* y,const double* z,const double* x2,const  double* y2,const double* z2,double* d, int size, int size2);
double distance(double x,double y,double z,double x1,double y1,double z1);
int argmin(const double* d,int s);
void zero(double* d);


void init_vec(vec & d,int N);

int main(int argc, char *argv[]) {


  int my_rank;
  int ierr;
  int comm_size;

  ierr = MPI_Init ( &argc, &argv );
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &comm_size );
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );


    int Npoints,Nbins,Npoints_total;

    Npoints_total=atoi(argv[1]);
    Npoints=(int) Npoints_total/comm_size;//atoi(argv[1])
    Nbins=atoi(argv[2]);


    double* x;
    double* y;
    double* z;

    double* x2;
    double* y2;
    double* z2;
    double* d;
    double* result;

    x= (double *) malloc(Npoints*sizeof(double));
    y= (double *) malloc(Npoints*sizeof(double));
    z= (double *) malloc(Npoints*sizeof(double));

    read_file("../Data/all.dat",x,y,z,my_rank,comm_size,Npoints);

    x2= (double *) malloc(Nbins*sizeof(double));
    y2= (double *) malloc(Nbins*sizeof(double));
    z2= (double *) malloc(Nbins*sizeof(double));
    d=  (double *)  malloc(Nbins*sizeof(double));
    
    //std::cout<<"start "<<my_rank<<' '<<x[0]<<' '<<y[0]<<' '<<z[0]<<' '<<d[0]<<std::endl;

    read_file("../Data/centers_6.dat",x2,y2,z2,0,1,Nbins);
    
    dis(x,y,z,x2,y2,z2,d,Npoints/comm_size,Nbins);

    result= (double *) malloc(Nbins*sizeof(double));
    
     ierr = MPI_Reduce(d,result,Nbins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (my_rank==0){
      std::ofstream outfile;
    outfile.open ("histogram.dat");
    outfile.precision(16); outfile.setf(std::ios::scientific);

        
    for (int j = 0; j<Nbins;j++){
    outfile<<x2[j]<<' '<<y2[j]<<' '<<z2[j]<<' '<<result[j]<<std::endl;
    }
        outfile.close();
    }
    
  MPI_Finalize();
  
    return 0;
}

void read_file(std::string filename,double * x, double * y,double * z,int id, int size, int size2){

    
    double xfile,yfile,zfile;
    
    //std::cout<<"Reading file "<<id<<' '<<size<<' '<<size2<<' '<<size2/size<<std::endl;
    std::stringstream ss;
    ss<<filename;
    std::fstream myfile(ss.str(), std::ios_base::in);
    int j=0;
    for (int i = 0; !myfile.eof() ; i++){
        if(i>=(id+1)*size2/size){
            myfile.close();
            break;
        }else if (i>id*size2/size){
            
        myfile >> xfile;
        myfile >> yfile;
        myfile >> zfile;

        x[j]=xfile;
        y[j]=yfile;
        z[j]=zfile;    
        j++;
        }else{
                myfile >> xfile;
                myfile >> yfile;
                myfile >> zfile;
        }
    }
}


void dis(const double* x ,const double* y ,const double* z,
         const double* x2,const double* y2,const double* z2,double* d, int size, int size2){

    double* dd;
    dd= (double *) malloc(size2*sizeof(double));   
    
    for (int i = 0; i<size;i++){    
        for (int j = 0; j<size2;j++){
            dd[j]=distance(x[i],y[i],z[i],x2[j],y2[j],z2[j]);        
        }   
    d[argmin(dd,size2)]+=1;  
    }

}

void zero(double* d){
    int size= (*(&d + 1) - d);
    for (int i = 0; i<size;i++){
        d[i]=0;  
    }
}

double distance(double x,double y,double z,double x1,double y1,double z1){
return sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1));
}


int argmin(const double* d,int size){
    double d_min=1.e10;
    int d_min_arg=-1;

    for (int i = 0; i<size;i++){
        if(d_min>d[i]){
            d_min=d[i];
            d_min_arg=i;
        }
    }
    return d_min_arg;
}
