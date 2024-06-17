try:
    from moobench.optbase import *
    from moobench.optlib_scipy import ScipyOptimizationAlgorithm
except ImportError as error:
    print('An exception occurred: {}'.format(error))
    pass

from grillage.grillage_model import Grillage

from grillage.grilage_beam_anan import generate_grillage_analysis,BeamPropertAnAn
from grillage.grillage_model import TBeamProperty,LBeamProperty

class CallbackPropertyGetSetConnector(BasicGetSetConnector):
    def __init__(self,prop_get_set):
        self._prop_get_set = prop_get_set
    @property
    def value(self):
        return self._prop_get_set

    @value.setter
    def value(self,value):
        self._prop_get_set = value

class CallbackGetPropRatioConnector(BasicGetConnector):
    def __init__(self,num,denum):
        self._num =num
        self._denum = denum
    @property
    def value(self):
        return self._num/self._denum


class GrillageBeam_AnMod(AnalysisExecutor):
    def __init__(self,grillage:Grillage):
        super().__init__()
        self.grill_model = grillage
        self.grill_anan = generate_grillage_analysis(grillage)


    def analyze(self):
        self.grill_anan.calculate_and_apply_elastic_reactions()

        return AnalysisResultType.OK

class GrillageBeamOptimization_OptProb(OptimizationProblem):
    def __init__(self,grillage:Grillage,name=''):
        if name == '':
            name = 'Grillage Beam Optimization'
        super().__init__(name)
        am = GrillageBeam_AnMod(grillage)
        for beam in am.grill_anan.beams.values():
            prop:BeamPropertAnAn = beam.prop
            if prop.prop is TBeamProperty or prop.prop is LBeamProperty:
                tprop:TBeamProperty = prop.prop
                self.add_design_variable('hw_'+str(beam.id),CallbackPropertyGetSetConnector(tprop.hw),5.0,25.0)
                self.add_constraint('g_bf_tf_'+str(beam.id),CallbackGetPropRatioConnector(tprop.bf,tprop.tf),25.0,ConstrType.LT)
                self.add_constraint('g_bf_tf_' + str(beam.id), CallbackGetPropRatioConnector(tprop.bf, tprop.tf), 4.0,
                                    ConstrType.GT)

            self.add_constraint(DesignConstraint('g_Sig_'+str(beam.id),CallbackGetConnector(beam.getSigmax_crit), prop.prop.mat, ConstrType.LT))

        self.add_objective('mass',CallbackGetConnector(am.grill_model.getGrillageMass))
        self.add_analysis_executor(am)

def prepare_and_run_optimization(grillage:Grillage):
    out_folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    isExist = os.path.exists(out_folder_path)
    if not isExist:
        os.makedirs(out_folder_path)

    op = GrillageBeamOptimization_OptProb(grillage)

    # SciPy algorithms
    if True:
        opt_ctrl = {}
        op.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000', 'SLSQP', opt_ctrl)
        if True:
            sol = op.optimize()
            op.print_output()
        else:
            sol = op.optimize_and_write(out_folder_path)