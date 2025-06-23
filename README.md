<h1 style="font-size:3.2em">🧠 Design of SPR Biosensors using Deep Attentive Network (DANet)</h1>


A physics-guided deep learning framework for predicting optimal metal and graphene layer thicknesses in multilayered Surface Plasmon Resonance (SPR) fiber optic biosensors — based on target sensitivity and Full Width at Half Maximum (FWHM) and other properties.


<h3>📌 Problem Statement </h3>
Designing an SPR biosensor isn’t just about choosing a metal — it involves selecting the right combination of materials and layer thicknesses (e.g., Silver, Gold, Copper, Aluminum etc..). Each metal’s thickness can vary (0–50 nm) and graphene layers from 0 to 10 — resulting in:
⚠️ Millions of potential combinations!

Traditionally, this design process is done via time-consuming trial-and-error using physics simulations. Each configuration must be manually tested using Maxwell-based modeling — a laborious and computationally expensive process.



<h3> 🔁 Our Approach: Inverse Design </h3>
Instead of forward simulation (input → output), we tackle the inverse problem:

We leverage the Inverse Transfer Matrix Method — rooted in Maxwell’s equations and multilayer optics — to simulate biosensor behavior across thousands of material combinations and refractive indices. Using this high-fidelity dataset, we train a deep attention-based neural network (DANet) to learn complex nonlinear mappings between biosensor performance (sensitivity, FWHM, resonance λ) and its structural design (thicknesses of silver, gold, copper, aluminum, and number of graphene layers).



<h3> 💡 Why It Matters </h3>


Surface Plasmon Resonance (SPR)-based biosensors are foundational in:

🧬 Early cancer and virus detection

🌊 Water contamination monitoring

🍽️ Food safety and toxin detection

🩺 Real-time medical diagnostics

🧪 Drug discovery & biomolecular analysis

Enabling smarter biosensor design saves time, resources, and cost — accelerating real-world applications in healthcare, environment, and biotechnology.



<h3>🔬 Dataset Generation</h3>

We generated a high-fidelity dataset (~80,000 samples) using:

📐 Inverse Transfer Matrix Method (physics-based modeling of multilayer optics)

🧮 Maxwell’s equations applied to SPR structures

<h3>⚙️ MATLAB simulation varying:</h3>

Metal layer thicknesses (Silver, Gold, Copper, Aluminum)

Number of graphene layers

Sensing medium refractive index

Fiber and sensor parameters


<h2>📊 Inputs & Outputs </h2>

🔢 Input Features:

FWHM_nm – Full Width at Half Maximum

Sensitivity_nm_per_RIU – Sensitivity in nm/RIU

NumA – Numerical Aperture

Sensing_length_mm – Length of sensing region

n_sens – Sensing medium refractive index

Resonance_lambda_nm – Resonance wavelength


🎯 Output Targets:

Silver_thickness_nm

Gold_thickness_nm

Copper_thickness_nm

Aluminum_thickness_nm

Graphene_layers

<h3>🧠 Model: Deep Attention Network (DANet) </h3>

We used DANet for multi-output regression due to its strength in capturing nonlinear interactions in tabular data.

📈 Features:
✅ Multi-output regression
✅ Built-in feature scaling and inverse transform
✅ GPU-accelerated with PyTorch
✅ Clean training pipeline using PyTorch Tabular
✅ Training vs validation vs test loss curves
✅ Evaluation via MAE per output variable


</h3>🛠️ Tech Stack </h3>

Python 3.8+

PyTorch + pytorch-tabular

DANet (Deep Attention Network)

Pandas, NumPy, Matplotlib

scikit-learn

omegaconf


<h3>📦 Installation </h3>

pip install torch torchvision torchaudio
pip install pytorch-tabular omegaconf pandas numpy matplotlib scikit-learn
<h3>🚀 Run It </h3>
python DANet.py


<h3>📈 Evaluation & Metrics </h3>

MAE (Mean Absolute Error) reported for each metal and graphene output

Training, validation, and test loss curves plotted

Test loss highlighted for performance comparison


<h3>🔧 Benefits of Our System </h3>

✅ Automates inverse biosensor design
✅ Saves significant cost and time
✅ Enables rapid prototyping and experimentation
✅ Bridges optical physics with deep learning
✅ Reduces manual trial-and-error from hours/days to milliseconds


<h3>📜 License </h3>

MIT License

<h3>🙌 Credits </h3>

📘 Simulation Inspired By:
“Theoretical Analysis of Sensitivity Enhancement by Graphene Usage in Optical Fiber Surface Plasmon Resonance Sensors”
— A. A. de Melo et al., IEEE Transactions on Instrumentation and Measurement, 2018

“Sensitivity evaluation of a multi-layered surface plasmon resonance-based fiber optic sensor”
— B. D. Gupta, A. K. Sharma, IIT Delhi, 2004

<h3>🤖 ML Implementation: </h3>
DANet via PyTorch Tabular
https://pytorch-tabular.readthedocs.io/
