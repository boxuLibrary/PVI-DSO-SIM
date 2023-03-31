#ifndef VIEWER_H
#define VIEWER_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pangolin/pangolin.h>

#include <param.h>

class viewer
{
private:
    

    void drawMap(const std::vector<Eigen::Vector3d>& map, int color = 1 , int pointsize = 5);                // draw mappoints
    void drawCamera(const Eigen::Matrix3d Rwc, const Eigen::Vector3d twc);    // draw 3d camera model
    void drawTraj(const std::vector<Eigen::Vector3d>& traj);
    void drawCoordinate();
    void drawFrameWall(double d,  double h);
    void drawFrameGround(double d,  double h);

    pangolin::OpenGlMatrix GetCurrentOpenGLCameraMatrix(const Eigen::Matrix3d Rwc, const Eigen::Vector3d twc);     // transform eigen matrix to opengl matrix

public:
    viewer(const Param& params);
    void guiRun(std::string camfile, std::string mptfile, std::string localmptfile);
    ~viewer();

    Param param;
};



#endif /* VIEWER_H */
