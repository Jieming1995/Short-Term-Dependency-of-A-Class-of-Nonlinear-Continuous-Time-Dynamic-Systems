# Short-Term-Dependency-of-A-Class-of-Nonlinear-Continuous-Time-Dynamic-Systems

Dynamic systems usually have long-term dependency. Recent work, however, shows that the long-term memory models like recurrent neural networks and short-term memory models like feed-forward neural networks have comparable performance in approximating the behavior of dynamic systems. While the pioneering researchers try to understand the short-term dependency in discrete time dynamic systems, this paper focuses on a class of continuous time dynamic systems characterized by ordinary differential equations. In this class of continuous time dynamic systems, only the variable can be observed, and its derivatives with respect to time cannot be measured directly. By analyzing the observability of the continuous time dynamic systems and the sampled systems, we show that the current output only depends on finite steps of history information when sampling fast. If the sampling is fast enough such that the Euler approximation can approximate the sampled system well, the current output only relies on the most recent $n$ steps of history information. Here, $n$ is the order of the ordinary differential equation. Later, we verified our results with the NARMAX method.

# In the simulation, we verify our main results using the NARMAX method in three examples. Here are descriptions of our codes.  
## 1.We use Simulink to generate the output sequences of each example. 

For example 1,2,3, the Simulink files are eg1.slx eg2.slx eg2.slx which are included in data folder. The output sequence of three examples with different sampling period and noise are put in the data folder too. For example, the eg3_001_00001.mat file contains the input and output sequence of examples 3 using sampling period 0.01s with noise variance 0.0001.  

## 2.The general NARMAX algorithm/function is created in the NARMAX.m.  

## 3. The SERR, FPE, and ERR values of each term are generated by NARMAX_eg1.m NARMAX_eg2.m NARMAX_eg3.m for three examples.  
## 4. The comparison of continuous time system and sampled system using different sampling period are given in files Outputs_eg1.m Outputs_eg2.m Outputs_eg3.m for three examples.  
