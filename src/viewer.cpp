#include <unistd.h>

#include "viewer.h"
#include "utilities.h"

viewer::viewer(const Param &param_) : param(param_) {
}

viewer::~viewer() {
}

pangolin::OpenGlMatrix viewer::GetCurrentOpenGLCameraMatrix(const Eigen::Matrix3d Rwc, const Eigen::Vector3d twc) {
    pangolin::OpenGlMatrix M;
    {
        M.m[0] = Rwc(0, 0);
        M.m[1] = Rwc(1, 0);
        M.m[2] = Rwc(2, 0);
        M.m[3] = 0.0;

        M.m[4] = Rwc(0, 1);
        M.m[5] = Rwc(1, 1);
        M.m[6] = Rwc(2, 1);
        M.m[7] = 0.0;

        M.m[8] = Rwc(0, 2);
        M.m[9] = Rwc(1, 2);
        M.m[10] = Rwc(2, 2);
        M.m[11] = 0.0;

        M.m[12] = twc[0];
        M.m[13] = twc[1];
        M.m[14] = twc[2];
        M.m[15] = 1.0;
    }
    return M;
}


void viewer::drawCamera(const Eigen::Matrix3d Rwc, const Eigen::Vector3d twc) {

    pangolin::OpenGlMatrix Twc = GetCurrentOpenGLCameraMatrix(Rwc, twc);
//    std::vector<GLdouble > Twc = {R(0,0), R(1,0), R(2,0),0,     R(0,1), R(1,1), R(2,1),0 ,     R(0,2), R(1,2), R(2,2),0,   t[0],t[1],t[2],1};

//    std::cout << R << std::endl << t <<std::endl;
    const float &w = 0.5;
    const float h = w * 0.75;
    const float z = w * 0.6;
    glPushMatrix();
    glMultMatrixd(Twc.m);

    glLineWidth(1);
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(w, h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(w, -h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(-w, -h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(-w, h, z);

    glVertex3f(w, h, z);
    glVertex3f(w, -h, z);

    glVertex3f(-w, h, z);
    glVertex3f(-w, -h, z);

    glVertex3f(-w, h, z);
    glVertex3f(w, h, z);

    glVertex3f(-w, -h, z);
    glVertex3f(w, -h, z);
    glEnd();

    glPopMatrix();

}

void viewer::drawTraj(const std::vector<Eigen::Vector3d> &traj) {
    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(0.0, 1.0, 0.0);
    for (size_t i = 0; i < traj.size(); i++) {
        glVertex3d(traj[i].x(), traj[i].y(), traj[i].z());
    }

    glEnd();

//    glPointSize(10);
//    glBegin(GL_POINTS);
//    glColor3f(0.0,1.0,0.0);
//    glVertex3d( traj[0].x(), traj[0].y(), traj[0].z());
//    glEnd();
}

void viewer::drawCoordinate() {

    //  设置大小
    glLineWidth(4);
    //  开始
    glBegin(GL_LINES);
    // 绿色
    // 设置颜色
    glColor3f(0, 1, 0);
    // 设置起点、终点坐标
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1, 0);
    //  蓝色
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1);
    // 红色
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(1, 0, 0);

    // 结束
    glEnd();

}

void viewer::drawFrameWall(double d, double h) {
    //  设置大小
    glLineWidth(4);
    //  开始
    glBegin(GL_LINES);

    Eigen::Vector3d vetex(-d, h, -d);
    Eigen::Vector3d vetex1(d, h, -d);
    Eigen::Vector3d vetex2(d, h, d);
    Eigen::Vector3d vetex3(-d, h, d);


    Eigen::AngleAxisd rot(-30. / 180. * M_PI, Eigen::Vector3d::UnitY());
    vetex = rot.toRotationMatrix() * vetex;
    vetex1 = rot.toRotationMatrix() * vetex1;
    vetex2 = rot.toRotationMatrix() * vetex2;
    vetex3 = rot.toRotationMatrix() * vetex3;



    // 红色1
    glColor3f(1, 0, 0);
    glVertex3f(vetex1.x(),  vetex1.z(), h);
    glVertex3f(vetex.x(),  vetex.z(), h);
    // 红色2
    glVertex3f(vetex1.x(), vetex1.z() , -h);
    glVertex3f(vetex.x(), vetex.z() , -h);
    // 红色3
    glVertex3f(vetex2.x(), vetex2.z() , h);
    glVertex3f(vetex3.x(), vetex3.z() , h);
    // 红色4
    glVertex3f(vetex2.x(), vetex2.z() , -h);
    glVertex3f(vetex3.x(), vetex3.z() , -h);
    // 红色5

    glVertex3f(vetex3.x(), vetex3.z(), h);
    glVertex3f(vetex.x(), vetex.z(), h);

    // 红色6
    glVertex3f(vetex3.x(), vetex3.z(), -h);
    glVertex3f(vetex.x(), vetex.z(), -h);
    // 红色7
    glVertex3f(vetex2.x(), vetex2.z(), h);
    glVertex3f(vetex1.x(), vetex1.z(), h);
    // 红色8
    glVertex3f(vetex2.x(), vetex2.z(), -h);
    glVertex3f(vetex1.x(), vetex1.z(), -h);


    glVertex3f(vetex.x(), vetex.z(), -h);
    glVertex3f(vetex.x(), vetex.z(), h);
    // 红色6
    glVertex3f(vetex2.x(), vetex2.z(), -h);
    glVertex3f(vetex2.x(), vetex2.z(), h);
    // 红色7
    glVertex3f(vetex1.x(), vetex1.z(), -h);
    glVertex3f(vetex1.x(), vetex1.z(), h);
    // 红色8
    glVertex3f(vetex3.x(), vetex3.z(), -h);
    glVertex3f(vetex3.x(), vetex3.z(), h);

    // 结束
    glEnd();

}


void viewer::drawFrameGround(double d, double h) {
    //  设置大小
    glLineWidth(4);
    //  开始
    glBegin(GL_LINES);

    // 红色1
    glColor3f(1, 0, 0);

    // 红色2
    glVertex3f(-d, d, h);
    glVertex3f(d, d, h);
    // 红色3

    glVertex3f(-d, -d, h);
    glVertex3f(d, -d, h);

    // 红色4
    glVertex3f(d, -d, h);
    glVertex3f(d, d, h);

    // 红色5
    glVertex3f(-d, -d, h);
    glVertex3f(-d, d, h);

    // 结束
    glEnd();

}

void viewer::drawMap(const std::vector<Eigen::Vector3d> &map, int color, int pointsize) {

    glPointSize(pointsize);
    glBegin(GL_POINTS);
    if (color == 1)
        glColor3f(0.0, 1.0, 0.0);
    else
        glColor3f(0.0, 0.0, 1.0);

    for (size_t i = 0; i < map.size(); i++) {
        glVertex3d(map[i].x(), map[i].y(), map[i].z());
    }

    glEnd();


}

void viewer::guiRun(std::string camfile, std::string mptfile, std::string localmptfile) {

    bool setting_render_display3D = true;

    int w = 640;
    int h = 480;
    const int UI_WIDTH = 180;
    bool running = true;

    pangolin::CreateWindowAndBind("Main", 2 * w, 2 * h);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    glEnable(GL_DEPTH_TEST);

    // 3D visualization
    pangolin::OpenGlRenderState Visualization3D_camera(
            pangolin::ProjectionMatrix(w, h, 400, 400, w / 2, h / 2, 0.1, 500),
            pangolin::ModelViewLookAt(0.0, 0.0, -10., 0, 0, 0, pangolin::AxisY)
//            pangolin::ModelViewLookAt(0, 0.0, 0, 0,0,0, pangolin::AxisZ)
    );

    pangolin::View &Visualization3D_display = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(UI_WIDTH), 1.0, -w / (float) h)
            .SetHandler(new pangolin::Handler3D(Visualization3D_camera));

    // parameter reconfigure gui
    pangolin::CreatePanel("ui").SetBounds(0.0, 1.0, 0.0, pangolin::Attach::Pix(UI_WIDTH));


    pangolin::Var<bool> settings_showCurrentCamera("ui.CurrCam", true, true);
    pangolin::Var<bool> settings_showTrajectory("ui.Trajectory", true, true);
    pangolin::Var<bool> settings_Follow("ui.Follow", false, true);
    pangolin::Var<bool> settings_show3D("ui.show3D", true, true);

    // load data to draw
    std::vector<MotionData> campose;
    LoadPose(camfile, campose);
    std::vector<Eigen::Vector3d> traj;
    for (size_t i = 0; i < campose.size(); i++) {
        traj.push_back(campose[i].twb);
    }
    std::vector<Eigen::Vector3d> allmappoints;
    LoadPoints(mptfile, allmappoints);

    // Default hooks for exiting (Esc) and fullscreen (tab).
    int i = 0;
    while (!pangolin::ShouldQuit() && running) {
        // Clear entire screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (setting_render_display3D) {

            Eigen::Matrix3d Rwc = campose.at(i).Rwb;
            Eigen::Vector3d twc = campose.at(i).twb;

            pangolin::OpenGlMatrix Twc = GetCurrentOpenGLCameraMatrix(Rwc, twc);
            // 在机器人后方放置一个摄像机，这个外参数是自己通过调试找到的，调试方法是打印矩阵 Visualization3D_camera.GetModelViewMatrix()
            // Eigen::Matrix4d Tvr_eigen;
            // Tvr_eigen << 0, -1, 0, -0.0,
            //                 0, 0, 1, -1.,
            //                 -1, 0, 0, -3.,
            //                 0, 0, 0, 1;
            // pangolin::OpenGlMatrix Tvr(Tvr_eigen);
            // //  Tvw = (Twr * Trv).inverse() = (Twr * (Tvr * Trr).Inverse() ).Inverse();
            // Tvw = (Twr * (Tvr * Trr).Inverse() ).Inverse();  // model view is Tvw, means world to camera view frame.
            if (settings_Follow.Get())
                Visualization3D_camera.Follow((Twc).Inverse());

            // Activate efficiently by object
            Visualization3D_display.Activate(Visualization3D_camera);

            std::stringstream ss;
            ss << localmptfile << i << ".txt";
            std::vector<Eigen::Vector3d> pobs;
            std::vector<Eigen::Vector3d> localmpts;
            LoadPointObs(ss.str(), localmpts, pobs);

            drawMap(allmappoints, 0, 5);
            // drawMap(localmpts, 1, 10);
            // drawCamera(Rwc, twc);
            drawCoordinate();

            if (param.constraint_type == 0) {
                drawFrameWall(7, 2);
            } else {
                drawFrameGround(7, 5);
            }

            drawTraj(traj);

        }

        setting_render_display3D = settings_show3D.Get();

        // Swap frames and Process Events
        pangolin::FinishFrame();

        usleep(50000);
        i++;
        if (i > campose.size() - 1) {
            i = 0;
        }

    }


    printf("QUIT Pangolin thread!\n");
    printf("I'll just kill the whole process.\nSo Long, and Thanks for All the Fish!\n");

    exit(1);
}