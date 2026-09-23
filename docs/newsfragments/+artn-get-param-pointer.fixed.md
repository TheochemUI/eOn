get_param and get_runparam now take a void** out pointer. pARTn allocates the value and writes that pointer through the argument, so a void* out parameter stored it at the wrong address.
