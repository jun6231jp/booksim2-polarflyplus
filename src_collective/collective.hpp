#ifndef _COLLECTIVE_HPP_
#define _COLLECTIVE_HPP_

//#define PFP_CYCLE_DEBUG // tarfficmanager
//#define HCUBE_CYCLE_DEBUG 
//#define FATTREE_CYCLE_DEBUG

//#define PFP_STEPUP_DEBUG
//#define HCUBE_STEPUP_DEBUG
//#define FATTREE_STEPUP_DEBUG
//#define HCUBE_ROUTING_DEBUG //routefunc
//#define TORUS_ROUTING_DEBUG //routefunc
//#define FATTREE_ROUTING_DEBUG //routefunc
//#define FATTREE_DEBUG
//#define SWITCH // for 1-tier fat-tree

#define PAIRWISE
//#define RING
//#define BYNARYTREE

#ifdef PAIRWISE
const int max_step = 20;
#endif

#ifdef RING
const int max_step = 8192;
#endif

#ifdef BYNARYTREE
const int max_step = 20;
#endif

//#define PAIRWISE_TRAFFIC_DEBUG // traffic
//#define RING_TRAFFIC_DEBUG // traffic

#endif

