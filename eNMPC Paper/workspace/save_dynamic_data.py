from pyomo.environ import (
    ComponentUID,
)
from pyomo.environ import value as pyo_val
import json


def initialize_record_comp_values(variables, init_t=None):
    record = dict()
    for var in variables:
        if var.is_reference():
            cuid = ComponentUID(var.referent)
        else:
            cuid = ComponentUID(var)

        key = str(cuid)
        if init_t is None:
            record[key] = []
        else:
            record[key] = [pyo_val(var[init_t])]
    return record

def record_comp_values(mod, tp, record):
    for key, valuelist in record.items():
        cuid = ComponentUID(key)
        var = mod.find_component(cuid)
        if tp is None:
            valuelist.append(pyo_val(var))
        else:
            valuelist.append(pyo_val(var[tp]))
    return record

def save_dict_in_json(dictData, filePath):
    with open(filePath, "w") as fwrite:
        json.dump(dictData, fwrite)

# path = join(folder, "state_record.json")
# save_dict_in_json(state_record, path)

# path = join(folder, "control_record.json")
# save_dict_in_json(control_record, path)

# path = join(folder, "economic_cost_record.json")
# save_dict_in_json(economic_cost_record, path)

# path = join(folder, "cur_LyaFun_record.json")
# save_dict_in_json(cur_LyaFun_record, path)

# path = join(folder, "pre_stageCost_record.json")
# save_dict_in_json(pre_stageCost_record, path)

# path = join(folder, "notOptimal_step_record.json")
# save_dict_in_json(notOptimal_step_record, path)



# path = join(folder, "var_initialize_iter41.json")
# with open(path, 'r') as fr:
#     var_info = json.load(fr)
# for comp in olnmpc.component_data_objects(Var):
#     comp.set_value(var_info[str(ComponentUID(comp))])

# path = join(folder, "Lya_info_iter41.json")
# with open(path, 'r') as fr:
#     Lya_info = json.load(fr)
# for key, val in Lya_info.items():
#     var = olnmpc.find_component(key)
#     var.set_value(val)

# path = join(folder, "curr_state_iter41.json")
# with open(path, 'r') as fr:
#     curr_states = json.load(fr)