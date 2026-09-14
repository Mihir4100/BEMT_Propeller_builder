# BEMT_Propeller_builder
A MATLAB-based propeller design and analysis toolkit utilizing Blade Element Momentum (BEM) theory. It calculates optimized sectional twist and chord distributions from airfoil polars and predicts thrust and torque performance.

How to use 
1. Extract .zip to get all the .csv files for S1223 airfoil
2. Make sure all the files are in same directory
3. open Blade_Designer.m and s1223_polars.m in matlab
4. Run S1223_polars.m once. It will create a S1223_polars.mat files which contains Aerodynamic performance of S1223 airfoil at diffrent reynolds number.
5. Now Run Blade_Designer.m as per required conditions which will output a table for different sections and also give predicted Thrust.
Note: Tip loss function can be toggled on and off 
