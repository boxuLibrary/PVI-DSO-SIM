# PVI-DSO-SIM

A basic implementation of coplanar constraints used in PVI-DSO [PVI-DSO: Leveraging Planar Regularities for Direct Sparse Visual-Inertial Odometry](https://arxiv.org/abs/2204.02635) in simulation environment.


<!-- CodeStructure -->
### CodeStructure
&emsp;&emsp;We release our basic code structure. In the code, we provide two scenarios (wall and ground) to test the different plane constraints, which can be seen intuitiely from the visual interface.
Meanwhile, the constraints of three different modes described in the paper are implemented, including: `Points`, which represents the point-based method; `PP (-L)` which represents 
fusing the co-planar constraints with the loosely coupled method; `PP (-T)` which represents fusing the co-planar constraints with tightly coupled method. It should be noted that
to improve the efficiency of optimization, we do not use automatic derivation or numerical derivation. Linearization of all the constraints is performed by <font color="Blue">analytical Jacobian</font>. 

### How to run
#### 1. Prerequisites
1.1 **Ubuntu** and **python**

* Ubuntu 16.04 or Ubuntu18.04

1.2. **Dependency**

* C++14 or C++17 Compiler
* Eigen 3.3.7
* OpenCV 3.4.9
* Boost 1.58.0
* Cere-solver 1.14.0: [Ceres Installation](http://ceres-solver.org/installation.html), remember to **sudo make install**.

#### 2. Build Project with Cmake
Clone the repository and compile the project:
```
git clone https://github.com/boxuLibrary/PVI-DSO-SIM.git
cd ~/PVI-DSO-SIM-master/
mkdir build
cd build
cmake ..
make -j4
```
#### 3. Run program
**Notice**: The executable file **data_gen** is in the bin directory, and you can run it by **./data_gen**

```
cd ../bin
./data_gen
```

<!-- CONTRIBUTING -->
### Related Papers

- **PVI-DSO: Leveraging Planar Regularities for Direct Sparse Visual-Inertial Odometry**.
```
@article{xu2022pvi,
  title={PVI-DSO: Leveraging Planar Regularities for Direct Sparse Visual-Inertial Odometry},
  author={Xu, Bo and Li, Xin and Li, JianCheng and Yuen, Chau and Dai, JiCheng and Gong, YiQun},
  journal={arXiv preprint arXiv:2204.02635},
  year={2022}
}
```
If you use the code for your academic research, please cite our related papers. And we use [vio_data_simulation ](https://github.com/HeYijia/vio_data_simulation.git) 
as our basic code and generate the simulation data.
<!-- LICENSE -->
### License

Distributed under the MIT License.


<!-- CONTACT -->
### Contact
We are still working on improving the code reliability. For any technical issues, please feel freely to contact us:

Bo Xu  - boxu1995@whu.edu.cn.

Xin Li - lixin.1997.lixin@gmail.com








