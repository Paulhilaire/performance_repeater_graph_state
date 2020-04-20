from fiber import Fiber
from tree_cluster import TreeCluster
from repeater_graph_state import RepeaterGraphState
from error_model import SingleQubitError
from rgs_network import RGSNetwork
from source import SinglePhotonSource
from detector import Detector
from rate import Rate

# Parameters of the protocol.
param = {
    "T_CZ": 1,  # CZ gate time
    "L": 50,  # Total distance
    "m": 14,  # Number of arms of an RGS (divided by two)
    "L_0": 0.19,  # Separation distance between source nodes
    "branches": [10, 5],  # branching parameter of the RGS
    "epsilon": 0.00002,  # single photon qubit error rate
    "eta_c": 0.98,  # collection efficiency of single photons
    "eta_d": 0.98,  # detection efficiency of single photons
}


# #############################################################
# Defining all the components of the network.

# Define a fiber
fiber = Fiber()

# Define a matter qubit single photon source.
source = SinglePhotonSource()

# Define a single photon detector.
detector = Detector()

# Define a error correction tree graph state
tree = TreeCluster()

# Define a repeater graph state
rgs = RepeaterGraphState()
rgs.tree = tree  # Attach the error correcting tree object to the repeater graph state

# Define a network of repeater graph states.
network = RGSNetwork()
network.fiber = fiber  # Attach the fiber object to this network.
network.RGS = rgs  # Attach the repeater graph state object to this network.
network.source = source  # Attach the source object to this network.
network.detector = detector  # Attach the detector object to this network.

# Define an error model (here a single qubit error model)
error_model = SingleQubitError()
error_model.network = network  # Attach the network object to this error model.

# Define the object in which the rate are calculated
rate = Rate()
rate.network = network  # Attach the network object to it.
rate.error = error_model  # Attach the error model object to it.


# #############################################################
# Place the protocol parameters.

# Define the branching parameter of the error correction tree graph state
tree.branches = param["branches"]

# Define the single photon collection efficiency
source.proba = param["eta_c"]

# Define the single photon detection efficiency
detector.det = param["eta_d"]

# Define the CZ gate time for the generation of the RGS
rgs.T_CZ = param["T_CZ"]
# Define the number of arms of the repeater graph state
rgs.m = param["m"]

# Define a total distance to this protocol
network.distance = param["L"] * fiber.attenuation_distance
# Define the separation distance between source nodes
network.L_0 = param["L_0"] * fiber.attenuation_distance

# Define the single qubit error.
error_model.epsilon = param["epsilon"]


# #############################################################
# Calculate the rates.
rate.get_rates()

# All the results are stored in this object. (examples)

# Probability of success of the protocol.
print(rate.proba)

# Rate of the protocol (probability / generation time of an RGS (rgs.T)) / number of matter qubits.
print(rate.rate_ss)

# Secret key rate per channel use   / per photon.
print(rate.skr_ch)
print(rate.skr_ph)
