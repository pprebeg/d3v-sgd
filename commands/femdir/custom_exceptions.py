class FeasibilityTestFailed(Exception):
    def __init__(self):
        super().__init__()

    def __str__(self):
        return "Grillage model does not pass mesh feasibility test! " \
               "Check if there are any plating zones between two adjacent" \
               " Primary Supporting Members with the same stiffener orientation" \
               " and different stiffener spacing."


class InvalidDesiredAspectRatio(Exception):
    def __init__(self, desired, maximum):
        super().__init__()
        self._desired = desired
        self._maximum = maximum
        self._message1 = "Desired plate element aspect ratio can not be greater" \
                         " than the maximum plate aspect ratio. Entered desired " \
                         "aspect ratio:"
        self._message2 = "Entered maximum aspect ratio:"

    def __str__(self):
        return f"{self._message1} {self._desired} {self._message2} {self._maximum}"


class NegativeRemainingDistance(Exception):
    def __init__(self, plate_id):
        super().__init__()
        self._plate_id = plate_id
        self._message1 = "Transition element dimension on plating zone"
        self._message2 = "is negative!, Primary supporting member flange " \
                         "overlaps a stiffener, check flange width values."

    def __str__(self):
        return f"{self._message1} {self._plate_id} {self._message2}"


class BaseMeshDimensionX(Exception):
    def __init__(self):
        super().__init__()
        self._message = " Can't return base mesh dimension x because the" \
                        " dictionary of base dimensions is empty! Calculate" \
                        " element base size mesh first."

    def __str__(self):
        return f"{self._message}"


class BaseMeshDimensionY(Exception):
    def __init__(self):
        super().__init__()
        self._message = " Can't return base mesh dimension y because the" \
                        " dictionary of base dimensions is empty! Calculate" \
                        " element base size mesh first."

    def __str__(self):
        return f"{self._message}"


class FeasibilityTestFailLong(Exception):
    def __init__(self, psm1_id, psm2_id, plate1_id, plate2_id, spacing1, spacing2):
        super().__init__()
        self._psm1_id = psm1_id
        self._psm2_id = psm2_id
        self._plate1_id = plate1_id
        self._plate2_id = plate2_id
        self._spacing1 = spacing1
        self._spacing2 = spacing2

    def message_string(self):
        line = str("Grillage model does not pass mesh feasibility test! ")
        line += str("Stiffener spacing of longitudinal stiffeners on plating "
                    "zones between adjacent longitudinal primary supporting members ")
        line += str(self._psm1_id)
        line += str(" and ")
        line += str(self._psm2_id)
        line += str(" do not match! Plating zone ")
        line += str(self._plate1_id)
        line += str(" has stiffener spacing of ")
        line += str(self._spacing1)
        line += str(", while plating zone ")
        line += str(self._plate2_id)
        line += str(" has stiffener spacing of ")
        line += str(self._spacing2)
        line += str(" m. Adjust these stiffener spacing values to match.")
        return line

    def __str__(self):
        return self.message_string()


class FeasibilityTestFailTran(Exception):
    def __init__(self, psm1_id, psm2_id, plate1_id, plate2_id, spacing1, spacing2):
        super().__init__()
        self._psm1_id = psm1_id
        self._psm2_id = psm2_id
        self._plate1_id = plate1_id
        self._plate2_id = plate2_id
        self._spacing1 = spacing1
        self._spacing2 = spacing2

    def message_string(self):
        line = str("Grillage model does not pass mesh feasibility test! ")
        line += str("Stiffener spacing of transverse stiffeners on plating "
                    "zones between adjacent transverse primary supporting members ")
        line += str(self._psm1_id)
        line += str(" and ")
        line += str(self._psm2_id)
        line += str(" do not match! Plating zone ")
        line += str(self._plate1_id)
        line += str(" has stiffener spacing of ")
        line += str(self._spacing1)
        line += str(", while plating zone ")
        line += str(self._plate2_id)
        line += str(" has stiffener spacing of ")
        line += str(self._spacing2)
        line += str(" m. Adjust these stiffener spacing values to match.")
        return line

    def __str__(self):
        return self.message_string()


class MeshV2FeasibilityFail(Exception):
    def __init__(self, p_dim_total, f_dim):
        super().__init__()
        self._p_dim_total = "{:.1f}".format(p_dim_total)
        self._f_dim = "{:.1f}".format(f_dim)

    def message_string(self):
        line = str("Mesh V2 feasibility test failed!")
        line += str(" Cannot generate FE mesh using Mesh Variant V2 with the "
                    "given mesh control parameters because element size is too small.")
        line += str(" Flange element with dimension ")
        line += str(self._f_dim)
        line += str("mm is greater than the sum of first two plating elements ")
        line += str("equal to ")
        line += str(self._p_dim_total)
        line += str("mm. ")
        line += str(" Increase the maximum plating or flange aspect ratios")
        line += str(" and try again, or select another meshing variant.")

        return line

    def __str__(self):
        return self.message_string()


class MeshV2FeasibilityFailGeneric(Exception):
    def __init__(self):
        super().__init__()

    @staticmethod
    def message_string():
        line = str("Mesh V2 feasibility test failed!")
        line += str(" Cannot generate FE mesh using Mesh Variant V2 with the "
                    "given mesh control parameters because element size is too small.")
        line += str(" Increase the maximum plating or flange aspect ratios")
        line += str(" and try again, or select another meshing variant.")

        return line

    def __str__(self):
        return self.message_string()
