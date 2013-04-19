import numpy

class MatrixStack:

    def __init__(self, A):
        """
        """
        self.data = A
        self.is_double = A.dtype in [numpy.float, numpy.float64, numpy.float128, numpy.float, numpy.float_]

    def add(self, M):
        if self.is_double:
          matrix_stack[double](self.data).iadd(matrix_view[double](M))
        else:
          matrix_stack[dcomplex](self.data).iadd(matrix_view[dcomplex](M))

    def sub(self, M):
        if self.is_double:
          matrix_stack[double](self.data).isub(matrix_view[double](M))
        else:
          matrix_stack[dcomplex](self.data).isub(matrix_view[dcomplex](M))

    def matmul_L_R(self, L, D, R):
        if self.is_double:
          matrix_stack[double](self.data).onsite_matmul_L_R(matrix_view[double](L), matrix_stack[double](D), matrix_view[double](R))
        else:
          matrix_stack[dcomplex](self.data).onsite_matmul_L_R(matrix_view[dcomplex](L), matrix_stack[dcomplex](D), matrix_view[dcomplex](R))

    def invert(self):
        if self.is_double:
          matrix_stack[double](self.data).invert()
        else:
          matrix_stack[dcomplex](self.data).invert()
