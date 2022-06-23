# GP_MRAC_Code
Code for testing a Model Reference Adaptive Control (MRAC) scheme to control the pitch of a plane under wing-rock dynamics, accounting for modeling errors by using Gaussian Process (GP) regression. Implementation based on [1].

## Disclaimer
This is an implementation of the method presented in [1] I presented for a Parameter Inference and State Estimation graduate class. I don't claim ownership over the methodology presented in [1].

## Files
* **GP_KL_MRAC_simple_test.m** : Code for testing a Model Reference Adaptive Control (MRAC) scheme to control the pitch of a plane under wing-rock dynamics, accounting for modeling errors by using Gaussian Process regression (hence, the algorithm is called GP-MRAC). In this case, the data in the training set will be discarded based on Kullback–Leibler (KL) divergence.
* **GP_OL_MRAC_simple_test.m** : Code for testing a Model Reference Adaptive Control (MRAC) scheme to control the pitch of a plane under wing-rock dynamics, accounting for modeling errors by using Gaussian Process regression (hence, the algorithm is called GP-MRAC). In this case, the oldest data point in the training set will be discarded.
* **RBFN_MRAC_simple_test.m** : Code for testing a Model Reference Adaptive Control (MRAC) scheme to control the pitch of a plane under wing-rock dynamics, accounting for modeling errors by using a Gaussian Radial Basis Function Network (RBFN).
* **MRAC_Ideal_Tuner.m** : Code for tuning Model Reference Adaptive Controller (MRAC) Gains (K), without taking into account modelling errors (ideal case). Thus, Wing-Rock Dynamics are represented by a double integrator plant
* **Reference_Model_Tuner.m** : Code for tuning reference model response for Model Reference Adaptive Control (MRAC) reference tracking.
* **Bayesian Nonparametric MRAC Using GPR_Final_Report.pdf** : 
* **Bayesian Nonparametric MRAC Using GPR_Presentation.ppt** : Power Point presentation I used to present the results from using GP-MRAC to control the pitch of a plane under wing-rock dynamic.
* **Project_Proposal.pdf** : Document I used to propose the GP-MRAC implementation as a final project for the Parameter Inference and State Estimation course to the instructor.

## Citing work
* **[1] G. Chowdhary, H. A. Kingravi, J. P. How, and P. A. Vela,**. "Bayesian nonparametric adaptive control using gaussian processes," IEEE Transactions on Neural Networks and Learning Systems, vol. 26, pp. 537-550, March 2015.
