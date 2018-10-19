---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: "Stochastic Reachability Toolbox"
---

SReachTools is an open-source MATLAB Toolbox for performing stochastic verification and reachability analysis.  

- Can you show me some examples of SReachTools working? 
    - We have cataloged a number of [examples ](https://unm-hscl.github.io/SReachTools/examples/) implemented using SReachTools. These examples are also available as part of the source code of SReachTools, see `examples/*.m`. 
- Where do I get the source code from? How do I install this? What are the dependencies?
    - See our [Github repository](https://github.com/unm-hscl/SReachTools). We do periodic releases regarding bug fixes which can be obtained from our [release page](https://github.com/unm-hscl/SReachTools/releases). Our [quick start guide](#quick-start-guide), described further down this page, walks through the installation process.
- How can I use this toolbox? What are the terms and conditions to follow to use SReachTools?
    - SReachTools is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See our [License](license/). 
- Can I to contribute to this toolbox?
    - Of course, we welcome pull requests. See [Contributing guidelines](contributing/). 
- Where do I ask questions or give feedback? 
    - Use our [Google groups](https://groups.google.com/d/forum/sreachtools) or the [Github issues](https://github.com/unm-hscl/SReachTools/issues) page.

The authors of this toolbox are [Abraham P.  Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.  Gleason](http://www.unm.edu/~gleasonj/). Please cite their [relevant papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using the toolbox. The authors are PhD advisees of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).

We have submitted a tool paper describing the features of SReachTools to the *22nd ACM International Conference on Hybrid Systems: Computation and Control summarizing the features of SReachTools*. A copy of this submission is [available in the repository](https://github.com/unm-hscl/SReachTools/raw/master/SReachTools.pdf).

{% include news.html %}

## What does SReachTools do?

SReachTools focuses on the following problem --- Construct **controllers** and characterize the **set of initial states** such that 
1. the controller satisfies the specified control bounds,
1. the stochastic system stays within a time-varying **target tube** with a probability above a given threshold? \\
For example, a typical **reach-avoid** constraint is to stay within a *safe set* to stay within the time horizon and reach a *target set* at the time horizon when starting from an initial state \\(\overline{x}\_0\\), as shown in the figure below.
<div class="desc-figure">
    <img src="{{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }}" alt="A cartoon depicting the stochastic reach-avoid problem"/>
    <!-- [A cartoon depicting the stochastic reach-avoid problem]({{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }})\\ -->
</div>
Here, we would like to pick the *green* controller over the *red* controller and compute the collection, the *orange set*, of all initial states such that the probability of success (reach-avoid) \\(\mathbb{P}\\) is above a given threshold \\(\theta\\).

This problem appears in a wide range of applications --- space applications ([spacecraft rendezvous and docking problem](./examples/FtCVXUnderapproxVerifyCWH.html)), transport (semi-autonomous/fully-autonomous cars and airplanes), biomedical applications (automated anesthesia delivery system), to name a few.

This toolbox provides MATLAB APIs to tackle this problem for Gaussian-perturbed linear time-invariant systems using [Fourier transforms](./examples/FtCVXUnderapproxVerifyCWH.html) [^1], [Lagrangian (set-operations)](./examples/lagrangianApproximations.html) [^2], and [dynamic programming](./examples/doubleIntegratorDynamicProgramming.html) [^3] [^4] methods.
We currently provide polytopic underapproximation and open-loop controller synthesis for this problem.
In future, we will provide extensions to linear time-varying systems, closed-loop controller synthesis, and non-Gaussian disturbances.

[^1]: A. P. Vinod and M. M. K. Oishi, "[Scalable Underapproximative Verification of Stochastic LTI Systems using Convexity and Compactness](https://doi.org/10.1145/3178126.3178148)", in Proceedings of Hybrid Systems: Computation and Control, 2018
[^2]: J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "[Underapproximation of Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian Methods](https://doi-org/10.1109/CDC.2017.8264291)," in Proceedings of the IEEE Conference on Decision and Control, 2017
[^3]: S. Summers and J. Lygeros, "[Verification of discrete time stochastic hybrid systems: A stochastic reach-avoid decision problem](https://doi.org/10.1016/j.automatica.2010.08.006)," Automatica, 2010.
[^4]: A. Abate, M. Prandini, J. Lygeros, S. Sastry, "[Probabilistic reachability and safety for controlled discrete time stochastic hybrid systems](https://doi.org/10.1016/j.automatica.2008.03.027)," Automatica, 2008.

## Quick start guide

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools.

1. MATLAB (>2017a)
    1. Toolboxes
        1. MATLAB's Statistics and Machine Learning Toolbox
        1. MATLAB's Global Optimization Toolbox (**Optional**)
1. MPT3 ([https://www.mpt3.org/](https://www.mpt3.org/))
    1. Copy the MATLAB script [install_mpt3.m](https://www.mpt3.org/Main/Installation?action=download&upname=install_mpt3.m) provided by MPT3 from the browser, and run it in MATLAB to automatically download MPT3 and its dependencies.
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/))
    1. Install the CVX (Standard bundle, including Gurobi and/or MOSEK)
    1. Installation instructions are given in [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).
1. (**Optional**) We recommend using Gurobi as the backend solver for the convex programs
   formulated by SReachTools. In practice, we find both CVX and MPT3 perform
   much better with Gurobi. See
   [http://www.gurobi.com/registration/download-reg](http://www.gurobi.com/registration/download-reg)
   for more details. Note that Gurobi offers free academic license.

### Installation

1. Install the necessary dependencies listed above
1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/unm-hscl/SReachTools/releases))
1. Change the MATLAB current working directory to where SReachTools was
   downloaded
1. Run `srtinit` in MATLAB to add the toolbox to the paths and ensure all
   must-have dependencies are properly installed.
   - You can add `cd <path_to_sreachtools_repo>;srtinit` to your MATLAB's
     `startup.m` to automatically have this done in future.
   - Additional steps (optional):
       - Run `srtinit -t` to run all the unit tests.
       - Run `srtinit -v` to visualize the steps the changes to the path and
         check for recommended dependencies.  
       - Run `srtinit -x` to remove functions of SReachTools from MATLAB's path
         after use.  

------