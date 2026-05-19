import matplotlib as mpl
import altair as alt

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']


def _arial_theme():
    return {'config': {'font': 'Arial'}}


alt.themes.register('arial', _arial_theme)
alt.themes.enable('arial')
