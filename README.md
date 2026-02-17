# SimScape-Electrical-Real-Time
Contains realistic power systems models made with SimScape Electrical, ready for parallel execution in real-time.

Since January 2026, the toolbox previously known as SimPowerSystems is not supported anymore by Mathworks.
The global objective of this repository is to provide user with equivalent power systems models in the native SimScape Electrical for Simulink.

SimScape Electrical was not designed by power systems engineers (as it is part of a much more general simulation tool) and lacks of some expertise in models, components and methods.
SimPowerSystems in contrast was made with real-time and power systems in mind by Hydro-Quebec.

Nevertheless, SimScape Electrical is a really powerful toolbox, only lacking some twists to make it dominent in the industry. The SimScape solver engine is a general multi-domain solver and is very similar to a nodal admittance solver like EMTP.

In this repository, we find:

- realistic power systems models, actually inspired by the SimPowerSystems in some case, often ready for real-time deployment in SpeedGoat using concurrent and parallel task processing: HVDC thyristor, HVDC-MMC, etc...
  
- some optimized model for real-time such as: half-bridge MMC, inverter switching functions, realistic decoupling elements such decouping distributed parameter lines.
  
- some very useful blocks for everyday modeling like configuable RLC block WITH MATCHING ICONS or unbalanced 3-phase RLC loads.

- Some model specifically designed for FPGA and HDL Coder (in /FPGA)

I don't have a SpeedGoat simulator so if you want to test the model and report the performance, that would be great. It is also possible to port them to other RTDS with some skills.

Folder description:

MATE: main Simulink library and associated files. Put in your MATLAB path to use the other models.
MATE/doc : contains various reading for real-time simulation aficionado.
MATE/slx/simetranss.slx: my dream library of models. For example, a unified configurable RLC series/parallel block with matching icon.
Note: if some models accidentally point to simetran.slx library (my first library), the corresponding model can be found in simetranss.slx.

Converter_MultiWinding_9LevelCascadedH-Bridge: a high power converter using shifted multi-winding transformer on the feeder to minimize grid harmonics and 9 level cascaded H-bridge 3-phase inverter. Optimimized for real-time implementation, with decoupled DC stages.

HVDC-MMC-200cellsPerArms: a full HVDC MMC bipolar link with distributed parameter lines with detailed optimized 200 cells per arms MMC model. Ready for concurent/parallel execution on SpeedGoat, using 4 different tasks.

HVDC_bipolar_Thyristor: a full HVDC bipolar link with distributed parameter lines, switched filter/reactive power compensation banks. Ready for concurent/parallel execution on SpeedGoat, using 4 different tasks.

NetworkDecouplingTechnique: The IEEE 13 node fedder SimScape model is used to demonstrate proper ways to decouple power system equation thougth the use of stublines and stublines transformer. The model and various techniques are explained in more details in MATE/doc/MATE_StublineHowTo.pdf

Stublines are Bergeron Line model with losses adjusted to have exactly one-time setp of propagation delay.
Stubline are commonly used in real-time simulators with power transformer. By replacing the secondary leakeage
inductance with a stubline, one obtains an accurate point of decoupling for the grid equations

DistributionGrids: Using the techniques of /NetworkDecouplingTechnique, distribution grid models are presented here.
List of models: IEEE 123 node test feeder model.

OversampledPowerConverters: Explains how to simulate high-frequency PWM converter with oversampling of IGBT/GTO/MOSFET gate signal. In the model, the 2-level inverter feed a 3-phase RL load with back-EMF. Sample time is 30 us and the PWM frequency is 3 kHz. However, the accuracy of the simulation is equivalent to a 2 us sample time, why? Because of oversampling and averaging of the gate pulses. In real-time simulator, this resolution can go as low as 10 nanoseconds with FPGA boards. The model is explained in more details in MATE/doc/MATE_OversamplingPowerConverters.pdf

Machines/  : contain many machine types (Synchronous, induction) designed with constant Jacobian (admittance), best for real-time simulation because the SimScape solver do not need to refactorize its equations for these node (assuming that the SimScape solver can optimally order its nodes).

FPGA/: In this folder will be found models suitable for FPGA simulation and HDL Coder. These models contains only Pejovic switches and Switching Function Inverters with full rectification and high impedance mode support.

List of models in FPGA/:

    - Basic Pejovic switch test model
    
    - Single-phase 2-level inverter with RL load, made with Pejovic switches and compensated for switching losses.
    
    - MATE_PMSMDynamicSwitches.slx : complete PMSM drive with compensated Pejovic switches (also called Dynamic Switch)
         This demo is derived from https://www.mathworks.com/help/hdlcoder/ug/generate-hdl-code-for-simscape-models-using-dynamic-switch-approximation.html
         Using compensated Dynamic Switches results in very accurate DC current simulation. By comparison, the original demo is completely innacurate for Idc (you need to add the current sensor at the battery). Accurate DC current are important in battery management systems. 
         Without compensation, the battery will see a greater energy outflow. This is the well-known problem spurious power losses of the Pejovic method. The power losses will increase with the PWM frequency. See reference section, IECON paper.
         
     - MATE_PermanentMagnetSynchronousMachine.slx: complete 9 kHZ PMSM Drive with AC-link rectifier. The inverter is a switching function with rectification capability and fixed admittance. The PMSM motor is of unprecedented formulation with the constant part of the mutual inductance directly into the SimScape solver. It also have a fixed admittance as seen by the main SimScape DAE solver. See /MATE/doc for the PMSM equations.
     The complete model is designed to be compatible with HDL Coder without modifications as all components have fixed admittance.

     -MATE_SwitchedReluctanceMotor64Drive.slx: complete SRM drive with switching functions. Extremely accurate and HDL Coder ready.

     -MATEInductionMachineDTC.slx: complete induction machine with direct torque control made with Pejovic compensated switches and Fixed referential IM, therefore well-suited for HDL coder.
     
     -MATE_3level_DriveOversampling.slx: a simple 3-level NPC drive with time-step averaging/high-impedance capable switching function. The model includes all common topology cases: w/wo neutral return, w/wo floating back-EMF source, battery disconnect, rectification, etc...

     More info in readme of the /FPGA section
     
         
Pejovic switches offer the advantage of having a constant admittance in both ON and OFF position. This constant admittance avoid the need for refactorisation of nodal equation (or state-space permutation matrix), something to avoid in FPGA because it is very time costly.
As I am not expert (yet) on MATLAB HDL coder,nor do I have FPGA boards, I would like to have feedback on the actual implementation of these models.

- Cite as Kyle Diller, SimScape Electrical models for real-time, 2026.
