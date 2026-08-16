echo "FIELD_COMMANDS=[lsort [info commands *field*]]"
echo "DATASET_COMMANDS=[lsort [info commands *dataset*]]"
echo "EXPORT_COMMANDS=[lsort [info commands *export*]]"
echo "GRADIENT_COMMANDS=[lsort [info commands *grad*]]"
echo "HELP_COMMANDS=[lsort [info commands *help*]]"
catch {help create_field} help_result
echo "CREATE_FIELD_HELP=$help_result"
catch {create_field -help} create_help_result
echo "CREATE_FIELD_DASH_HELP=$create_help_result"
foreach command {list_fields export_variables calculate_field_value} {
    catch {help $command} command_help
    echo "HELP_${command}=$command_help"
}
set dataset [load_file lin_state_0020_des.tdr -name Eq231Probe]
echo "DATASET=$dataset"
catch {list_fields -dataset $dataset} fields_result
echo "FIELDS=$fields_result"
foreach expression {
    {grad(<ElectrostaticPotential>)}
    {gradient(<ElectrostaticPotential>)}
    {ddx(<ElectrostaticPotential>)}
    {diff(<ElectrostaticPotential>,X)}
} {
    catch {create_field -dataset $dataset -name GradientProbe -function $expression} result
    echo "EXPRESSION=$expression RESULT=$result"
}
exit
