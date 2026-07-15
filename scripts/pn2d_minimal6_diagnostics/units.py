_FACTORS={("cm^-3","m^-3"):1e6,("m^-3","cm^-3"):1e-6,("V/cm","V/m"):1e2,("V/m","V/cm"):1e-2,("cm/s","m/s"):1e-2,("m/s","cm/s"):1e2,("A/cm^2","A/m^2"):1e4,("A/m^2","A/cm^2"):1e-4,("cm^2/(V s)","m^2/(V s)"):1e-4,("m^2/(V s)","cm^2/(V s)"):1e4,("cm^-1","m^-1"):1e2,("m^-1","cm^-1"):1e-2,("cm^-3*s^-1","m^-3*s^-1"):1e6,("m^-3*s^-1","cm^-3*s^-1"):1e-6}
def convert_value(value:float,source_unit:str,target_unit:str)->float:
    if source_unit==target_unit:return value
    try:return value*_FACTORS[(source_unit,target_unit)]
    except KeyError as error:raise ValueError(f"unsupported conversion: {source_unit} -> {target_unit}") from error
