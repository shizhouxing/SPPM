#ifndef CONFIG_H
#define CONFIG_H
using namespace std;

namespace Config {
    const int max_tracing_depth = 10;
    const int input_buffer_size = 1024;
    const int checkpoint_interval = 100;    
    const double alpha = 0.7;
    const double gamma = 0.5;
    const double initial_radius = 1e-5;
    const double epsilon = 1e-6;
}

#endif
