import matplotlib.pyplot as plt
from pyomo.repn import generate_standard_repn

def plot_figure(x_data,
                y_data,
                title = None,
                xlabel = None,
                ylabel = None,
                linewidth = 1.5,
                xrange = None,
                yrange = None,
                fontsize = None,
                figsize = None,
                grid = True,
                legend_loc = None,
                border = None,
                savename = None,):

    if figsize is not None:
        plt.figure(num = title, figsize = figsize)
    else:
        plt.figure(num = title)

    if fontsize is not None:
        plt.rcParams.update({'font.size':fontsize})

    numyline = len(y_data)
    if type(y_data) is dict:
        for key, val in y_data.items():
            plt.plot(x_data, val, linewidth = linewidth, label = key)
    else: # y_data is nested lists
        for ele in y_data:
            plt.plot(x_data, ele, linewidth = linewidth)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if xrange:
        plt.xlim(xrange)
    if yrange:
        plt.ylim(yrange)

    if grid:
        plt.grid()

    if type(y_data) is dict:
        if legend_loc is not None:
            plt.legend(loc = legend_loc)
        else:
            plt.legend()

    if border is not None:
        plt.subplots_adjust(left = border["left"],
                            right = border["right"],
                            bottom =  border["bottom"],
                            top = border["top"],
                            wspace = border["wspace"],
                            hspace = border["hspace"])

    if savename:
        plt.savefig(savename)

    plt.show()


def organize_pyomo_expr(expr):
    repn_form = generate_standard_repn(expr)
    return repn_form.to_expression()

