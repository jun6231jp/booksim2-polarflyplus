BookSim Interconnection Network Simulator
=========================================

BookSim is a cycle-accurate interconnection network simulator.
Originally developed for and introduced with the [Principles and Practices of Interconnection Networks](http://cva.stanford.edu/books/ppin/) book, its functionality has since been continuously extended.
The current major release, BookSim 2.0, supports a wide range of topologies such as mesh, torus and flattened butterfly networks, provides diverse routing algorithms and includes numerous options for customizing the network's router microarchitecture.

---

# Random traffic pattern simulation
Use files in src directory.

## 6D x F7 PolarFly+ with uniform traffic 
Enable the following line in polarfly_tables.hpp and rebuild.
``` polarfly_tables.hpp
//#define USE_TABLE_1x1
//#define USE_TABLE_7x3
//#define USE_TABLE_13x4
//#define USE_TABLE_31x6
#define USE_TABLE_57x8
```

Set the config file as follows:
``` config.txt
topology = polarflyplus;
k = 6; //hypercube
n = 8; //polarfly:F7=8, F5=6, F3=4, F2=3
routing_function = dim_order;
num_vcs = 6;
traffic = uniform;
```

## 0D x F5 PolarFly+ (PolarFly) with uniform traffic
Enable the following line in polarfly_tables.hpp and rebuild.
``` polarfly_tables.hpp
//#define USE_TABLE_1x1
//#define USE_TABLE_7x3
//#define USE_TABLE_13x4
#define USE_TABLE_31x6
//#define USE_TABLE_57x8  
```

Set the config file as follows:
``` config.txt
topology = polarflyplus;
k = 0; //hypercube
n = 6; //polarfly:F7=8, F5=6, F3=4, F2=3
routing_function = dim_order;
num_vcs = 6;
traffic = uniform;
```

# Collective communication pattern simulation
Use files in src_collective directory.

In the config file, set the packet size and injection rate as follows, ensuring their product equals 1:
``` config.txt
packet_size = 1;//16384; //8192; //4096; //2048; //1024;//512;//256;//128;//64;//32;//16;//8;
injection_rate = 1.0;//0.00006103515625; //0.0001220703125; //0.000244140625; //0.00048828125;//0.0009765625;//0.001953125;//0.00390625;//0.0078125;//0.015625;//0.03125;//0.0625;//0.125;
```

## 3D x F3 PolarFly+ with single-NIC pairwise exchange algorithm 
Enable the follwing line in collective.hpp.
``` collective.hpp
#define PFP_CYCLE_DEBUG
#define PAIRWISE
//#define RING
```

Enable the following line in polarfly_tables.hpp and rebuild.
``` polarfly_tables.hpp
//#define USE_TABLE_1x1
//#define USE_TABLE_7x3
#define USE_TABLE_13x4
//#define USE_TABLE_31x6
//#define USE_TABLE_57x8  
```

Set the config file as follows:
``` config.txt
topology = polarflyplus;
hypercubeport = 3; //hypercube port
polarflyport = 4; //polarfly port F7:8, F5:6, F3:4, F2:3
nic = 1;
num_vcs = 6;
traffic = pairwise;
routing_function = dim_order;
```

## 6D hypercube with 6-NIC pairwise exchange algorithm 
Enable the follwing line in collective.hpp.
``` collective.hpp
#define PFP_CYCLE_DEBUG
#define PAIRWISE
//#define RING
```
Enable the following line in polarfly_tables.hpp and rebuild.
```polarfly_tables.hpp
#define USE_TABLE_1x1
//#define USE_TABLE_7x3
//#define USE_TABLE_13x4
//#define USE_TABLE_31x6
//#define USE_TABLE_57x8
```

Set the config file as follows:
``` config.txt
topology = polarflyplus;
hypercubeport = 6; //hypercube port
polarflyport = 0; //polarfly port F7:8, F5:6, F3:4, F2:3
nic = 6;
num_vcs = 6;
traffic = pairwise;
routing_function = dim_order;
```

## star topology (single-tier 32-port switch) with 1-NIC pairwise exchange algorithm 
Enable the follwing line in collective.hpp and rebuild.
``` collective.hpp
#define FATTREE_CYCLE_DEBUG
#define PAIRWISE
//#define RING
#define SWITCH 
```
The "SWITCH" should only be enabled in this simulation. 

Set the config file as follows:
``` config.txt
topology = fattree;
k = 32; //ary
n = 1; //dim
nic = 1;
traffic = pairwise;
routing_function = nca;
```

## 0D x F2 PolarFly+ (PolarFly) with 2-NIC ring algorithm 
Enable the follwing line in collective.hpp.
``` collective.hpp
#define PFP_CYCLE_DEBUG
//#define PAIRWISE
#define RING
```
Enable the following line in polarfly_tables.hpp and rebuild.
```polarfly_tables.hpp
//#define USE_TABLE_1x1
#define USE_TABLE_7x3
//#define USE_TABLE_13x4
//#define USE_TABLE_31x6
//#define USE_TABLE_57x8
```

Set the config file as follows:
``` config.txt
topology = polarflyplus;
hypercubeport = 0; //hypercube port
polarflyport = 3; //polarfly port F7:8, F5:6, F3:4, F2:3
nic = 2;
num_vcs = 6;
traffic = ring;
routing_function = dim_order;
```

# Options

If you enable the following comments in polarfly_tables.hpp as needed, and debug messages are available. 

``` polarfly_tables.hpp
//#define PFP_DEBUG // trafficmanager
//#define DEBUG_DRAIN // trafficmanager
//#define PFP_MAP_DEBUG // trafficmanager
//#define PFP_NET_DEBUG // network
//#define PFP_FAULT_DEBUG // network
//#define PFP_ROUTING_DEBUG // routefunc
//#define PFP_ROUTER_DEBUG // router
```
Fault avoiding routing is also available.
In the config file, set the following parameter.

``` config.txt
link_failures = 1;
fail_seed=time;
```
The "link_failures" define the number of failure nodes in the network.


