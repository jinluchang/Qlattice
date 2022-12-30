from libcpp.string cimport string as std_string
from libcpp cimport bool

cdef extern from "qlat-utils/timer.h" namespace "qlat":

    int get_id_node()
    int get_num_node()

    long& verbose_level()

    double get_time()
    double& get_start_time()
    double& get_actual_start_time()

    double get_total_time()
    double get_actual_total_time()

    cdef cppclass Timer:
        long long flops
        Timer()
        Timer(const std_string& fname)
        void start()
        void start(bool is_verbose)
        void stop()
        void stop(bool is_verbose)
        @staticmethod
        void display(const std_string& tag)
        @staticmethod
        void autodisplay()
        @staticmethod
        void display_stack()
        @staticmethod
        void display_stack_always()
        @staticmethod
        void reset(long max_call_times_for_always_show_info)
        @staticmethod
        void fork(long max_call_times_for_always_show_info)
        @staticmethod
        void merge()
