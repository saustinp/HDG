function [added, multed] = mathOpsIntegrated(in1, in2)
    %#codegen
    % for code generation, preinitialize the output variable
    % data type, size, and complexity 
    added = 0;
    % generate an include in the C code
    coder.cinclude('flux.h');
    % evaluate the C function
    added = coder.ceval('flux', in1, in2); 
    multed = in1*in2; 
    end

    