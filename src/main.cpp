#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

#include <opencv2/opencv.hpp>

#include "grid.hpp"
#include "solveur.hpp"
#include "scheme.hpp"
#include "physics.hpp"

int lenght = 200;

void profile_x(Grid<U> &g, int row)
{
        int h = 200;
        cv::Mat img = cv::Mat(h,lenght,CV_8UC3, cv::Scalar(0,0,0));
        
        double scale = 100.0;
        
        for(int i=0; i<lenght; i++)
        {
                double val = scale*g.at(i,row)[0] +10;
                if(h-val>0 && h-val<h){
                img.at<cv::Vec3b>(h-val,i)[0] = 255;
                img.at<cv::Vec3b>(h-val,i)[1] = 255;
                img.at<cv::Vec3b>(h-val,i)[2] = 255;
                }
        }
        cv::imshow("rhox",img);
        cv::waitKey(10);
}

void profile_y(Grid<U> &g, int col)
{
        int h = 200;
        cv::Mat img = cv::Mat(h,lenght,CV_8UC3, cv::Scalar(0,0,0));
        
        double scale = 100.0;
        
        for(int i=0; i<lenght; i++)
        {
                double val = scale*g.at(col,i)[0] +10;
                if(h-val>0 && h-val<h){
                img.at<cv::Vec3b>(h-val,i)[0] = 255;
                img.at<cv::Vec3b>(h-val,i)[1] = 255;
                img.at<cv::Vec3b>(h-val,i)[2] = 255;
                }
        }
        cv::imshow("rhoy",img);
        cv::waitKey(10);
}

void draw(Grid<U> &g, int k, int save_img, cv::VideoWriter video)
{
        int h = g.Ny();
        int w = g.Nx();
        cv::Mat img = cv::Mat(h,w,CV_8UC3, cv::Scalar(0,0,0));
        
        double scale = 250;
        
        #pragma omp parallel for
        for(int i=0; i<h; i++)
        {
        for(int j=0; j<w; j++)
        {
                unsigned char val = std::min(255.,scale*g.at(j,i)[k]);
                
                img.at<cv::Vec3b>(i,j)[0] = val;
                img.at<cv::Vec3b>(i,j)[1] = 0;
                img.at<cv::Vec3b>(i,j)[2] = 0;
                
        }
        }
        cv::imshow("rho",img);
        cv::waitKey(10);
        
        if(save_img)
          video.write(img);
}

int main()
{	
	// Setup simulation with
	// - physical size x
	// - physical size y
	// - grid size x
	// - grid size y
	Simulation sim(lenght/10,lenght/10,lenght,lenght, "HLLC"); 
	
	// Init grid with two values
	// (rho, rho u, rho v, e)
	U u1(1.0,0.0,0.0,2500.0);
	U u2(0.125,0.0,0.0,250.);
	// Set central circular overdensity
	sim.init_circle(u1,u2,40);
	// Set boundaries
	sim.boundaries_extend();
	
	// Save output as a video
	cv::VideoWriter video("../output/outcpp.avi",cv::VideoWriter::fourcc('M','J','P','G'),30, cv::Size(lenght,lenght),true);
	
	int n=0;
	while(n<500)
	{
	   std::cout<<"Iteration "<<n++<<std::endl;
	   
	   sim.EDP2D(0.0005);	   
	   
	   Grid<U> g = sim.getGrid();
	   profile_x(g,lenght/2);
	   profile_y(g,lenght/2);
	   draw(g,0,0,video);
	}
	
	video.release();
	cv::destroyAllWindows();
	
	return  0;
}
