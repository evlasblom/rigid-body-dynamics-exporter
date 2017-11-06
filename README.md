# Rigid Body Dynamics Exporter
_An Approach using Differential Geometry and Lagrangian Dynamics_

## Introduction
The Rigid Body Dynamics Exporter (RBDE) can be used to export rigid body dynamics of an articulated body to a desired programming language. Currently, only Matlab is supported, but it can be easily extended.

Although many fast algorithms exist for computing the rigid body dynamics of robotic manipulators, this implementation differs in several ways. First of all, techniques from differential geometry are used to express the positions and motions in Euclidean space, based on [1][2]. This results in a much nicer set of equations that defines the linear and rotational dynamics as one, instead of separated in two sets of equations as is commonly done [3]. A similar approach exist, using so-called spatial vectors [4][5]. Though, in my understanding, the approach using differential geometry is much more elegant, powerful and expressive. Secondly, the dynamics are casted into a Lagrangian form. As such, this implementation does not aim for speed. Instead, the collection of functions primarily exists to be able to quickly build exporters for any kinematic or dynamic quantity, which can then be used for any kind of robotic manipulator. This makes it extremely useful for rapid prototyping of any kind of model-based control technique in robotics.

## Contents
This repository contains first and foremost a collection of functions for kinematics and screw theory, to accompany the book _A Mathematical Introduction to Robotic Manipulation_ by Murray et al. [1]. It contains most of the important functions, as well as several functions to describe connectivity from [4][5]. In additional, several exporters of kinematic and dynamic quantities are provided.

A number of very useful examples is included. This includes numerical integration of the dynamic models with and without constraints, nonlinear state estmation, simple balance control of a humanoid, performing torque control via Matlab to V-REP, etcetera. See the [examples](examples/README.md) for more information.

## Contributing
The Rigid Body Dynamics Exporter was built during and after my graduation and used in several research projects. Development is currently on hold, but feel free to contact me.

## Bibliography
[1]: Stramigioli, Stefano, and Herman Bruyninckx. Geometry and screw theory for robotics. Tutorial during ICRA 2001 (2001).

[2]: Richard M Murray, Zexiang Li, S Shankar Sastry, and S Shankara Sastry. A Mathematical Introduction to Robotic Manipulation. CRC press, 1994.

[3]: Spong, Mark W., Seth Hutchinson, and Mathukumalli Vidyasagar. Robot modeling and control. Vol. 3. New York: Wiley, 2006.

[4]: Featherstone, Roy. A Beginner's Guide to 6-D Vectors (Part 1). Robotics & Automation Magazine, IEEE 17.3 (2010): 83-94.

[5] Featherstone, Roy. A Beginner's Guide to 6-D Vectors (Part 2)[Tutorial]. Robotics & Automation Magazine, IEEE 17.4 (2010): 88-99.
