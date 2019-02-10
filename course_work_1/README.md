# Coursework 1 F&Q list


1. loading PLY: you can also use igl::readPLY and igl::readPLY and writePLY.
example of loading a .off file in the previous tutorials
[practices_flann/src/main.cpp#L31](https://github.com/smartgeometry-ucl/compM080-compGV18-2019/blob/master/tutorial_002/cpp/practices_flann/src/main.cpp#L31)

2. sub-sampling step in ICP: we suggest you start with non-subsampling step (use all input points if your machine allows). After your implementation looks work, you can try sub-sampling and answer item `4. Instead of directly aligning...` 

3. rejection step in ICP: a simple one is to check normal vectors' angle

4. how to do question 1 `In class, we discussed the solution to the optimization to minimize alignment error... Refine the derivation done in class to solve the optimization`: 
   You need to work out formulation in page 122 ( 02_Registration_term2_19.pdf/p122:Shape Matching Rotation) (with weight w_i this time.)

## Results reference

**bun045.ply TO bun000.ply**  
no noise/no subsampling/no rejection/point2point  
![Image](/course_work_1/imgs/45_to_00.gif)
![Image](/course_work_1/imgs/45_to_00_fin.jpg) 
![Image](/course_work_1/imgs/45_to_00_err.jpg)  


**bun090.ply TO bun000.ply**  
no noise/no subsampling/no rejection/point2point  
![Image](/course_work_1/imgs/90_to_00.gif)
![Image](/course_work_1/imgs/90_to_00_fin.jpg) 
![Image](/course_work_1/imgs/90_to_00_err.jpg)
