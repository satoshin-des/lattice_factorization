import ctypes, math, sys, os, platform

pf = platform.system()

def lat_fact(N: int, input_bit: int = 1, print_info: int = 1) -> list:
    """Integer factorization method using lattices.

    Args:
        N (int): Bit-size of RSA-type integer or RSA-type integer(Even if N is not RSA-type, this algorithm can run normally).
        input_bit (int, optional): N is bit-size(1) or integer(0). Defaults to 1.
        print_info (int, optional): Print progress information(1) or not(0). Defaults to 1.

    Returns:
        list: List of two (prime) factors of integer.
    """
    if (not input_bit) and math.log2(N) < 27:
        print("Input integer is too small.")
        sys.exit(0)
    else:
        if pf == 'Linux':
            LatFact = ctypes.cdll.LoadLibrary("./libfact.so")
        elif pf == 'Windows':
            os.add_dll_directory(os.getcwd())
            LatFact = ctypes.cdll.LoadLibrary('lat_fact.dll')
        LatFact.lattice_factorization.restype = ctypes.POINTER(ctypes.c_longlong)
        LatFact.lattice_factorization.argtypes = ctypes.c_int, ctypes.c_char_p, ctypes.c_int

        a = LatFact.lattice_factorization(input_bit, str(N).encode('utf-8'), print_info)
        return int(a[0]), int(a[1])
