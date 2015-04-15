# code for "Feedback Temperature Dependence and Equilibrium Climate Sensitivity"
# Jonah Bloch-Johnson, Raymond Pierrehumbert, Dorian Abbot
# this code is written in Julia (http://julialang.org/)
using PyPlot, Roots, DataFrames, Distributions

# ----------------
# settings for PyPlot (Julia interface with matplotlib)
# ----------------

rc("font", size=20)
rc("xtick.major", pad=8)
rc("ytick.major", pad=8)
rc("font", weight="light")
rc("axes", linewidth=1)
rc("mathtext", fontset="stixsans")
# Helvetica as matplotlib font: http://blog.olgabotvinnik.com/blog/2012/11/15/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/

# ----------------
# helper functions
# ----------------

# warming for quadratic model
function ΔT(F,λ,a)
    if a == 0
        -F / λ
    else
        dtrx = λ^2 - 4a * F
        dtrx >= 0. ? (-λ - sqrt(dtrx))/2a : Inf
    end
end

# warming for cubic model
function ΔT(F,λ,a,b)
    try
        sort(fzeros(ΔT -> F + λ * ΔT + a * ΔT^2 + b * ΔT^3, 0, 500))[1]
    catch
        Inf
    end
end

# ----------------
# Figure 1
# ----------------
figure(1,(10,8))

subplot_label_x = .025
subplot_label_y = .91
legend_font_size = 12

subplot(221)
    # helper functions/values
    λ = -.88
    N(T,a) = 3.71 + λ * (T - 287) + a * (T - 287)^2
    Ts = 285:.2:297

    # legend label
    plot([0,0],[1,1],label="\$\\lambda = $λ\$ \$W/m^2/K\$",linestyle="none")
    # curves of N(T)
    plot(Ts,[N(T,    0) for T=Ts],"k--",label="linear, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,0),1))K\$",linewidth=2)
    plot(Ts,[N(T,-.035) for T=Ts],"b",linewidth=2,label="\$a_C=-.035\$, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,-.035),1))K\$")
    plot(Ts,[N(T, .03) for T=Ts],"g",linewidth=2,label="\$a_M=.03\$, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,.03),1))K\$")
    plot(Ts,[N(T, .058) for T=Ts],"r",linewidth=2,label="\$a_H=.058\$, \$\\Delta T_{2x} = \$ ?")

    # x axis
    plot(Ts,[0 for T=Ts],"k",linewidth=2)

    # forcing arrow
    arrow(287,.3,0,3,shape="full",facecolor="k",head_width=.2,head_length=.2)
    text(285.7,1.5,"\$F_{2x}\$")

    # warming arrow
    arrow(287.3,-.2,3,0,shape="full",edgecolor="k",facecolor="k",head_width=.2,head_length=.2)
    text(288.1,-.9,"\$\\Delta T_{2x}\$")

    # ticks
    plot([287,287],[-.2,.2],"k",linewidth=2)
    text(286.4,-1,"\$T_0\$")
    plot([287+ΔT(3.71,λ,-.035),287+ΔT(3.71,λ,-.035)],[-.2,.2],"k",linewidth=2)
    plot([287+3.71/-λ,287+3.71/-λ],                      [-.2,.2],"k",linewidth=2)
    plot([287+ΔT(3.71,λ,.03),287+ΔT(3.71,λ,.03)],  [-.2,.2],"k",linewidth=2)

    # lambda
    plot([286,288],[3.71-λ,3.71+λ],linewidth=4,color="magenta")
    text(287.3,3.9,"\$\\lambda\$",color="magenta")

    # labels and axes
    annotate("a",xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    xlabel("\$T\$ (\$K\$)")
    xlim(Ts[1],Ts[end])
    xticks(285:2:297,285:2:297)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    ylim(-2,8)
    legend(loc="upper right",prop={"size"=>legend_font_size})

subplot(222)
    # helper functions/values
    # data for first four doublings from Meraner et al (2013)
    Ns = [0,-4.26,-9.09,-14.24,-19.18] + 27
    Ts = [0,2.79,6.7,12.46,22.68] + 287.72

    # data for 32xCO2 from T. Mauritsen, personal communication
    transient_Ns = [21.751, 16.974, 13.861, 11.937, 10.511, 9.4495, 8.3018, 8.7235, 7.0251, 6.9253, 5.6998]
    transient_Ts = [291.50, 295.54, 298.69, 301.17, 303.23, 305.10, 306.65, 308.19, 309.47, 310.59, 311.67]

    T_grid = -5:.1:25

    # legend label
    plot([1000,2000],[1000,1000],linestyle="none",label="ECHAM6")
    # curves/plots of N(T)
    plot(Ts,Ns,"bo",label="\$\\Delta T\$ vs. \$F_{32x} - F\$")
    plot(T_grid+287.72,[27 - 1.52 * T for T=T_grid],"k--",linewidth=2,label="linear (\$\\lambda_M = -1.52\$)")
    plot(T_grid+287.72,[27 - 1.52 * T + .03 * T^2 for T=T_grid],"r-",linewidth=2,label="quad (\$a_M = .03\$)")
    plot(transient_Ts,transient_Ns,"kx",label="\$32 \\times pCO_2\$ run")

    # x axis
    plot(T_grid+287.72,[0 for T=T_grid],"k-",linewidth=2)

    # forcing arrow
    arrow(287.72,2,0,23,shape="full",facecolor="k",head_width=.5,head_length=.9)
    text(283.821,12.0155,"\$F_{32x}\$")

    # ticks
    plot([287.72,287.72],[-1,1],"k-",linewidth=2)
    text(286.583,-4.34444,"\$T_0\$")

    # labels and axes
    annotate("b",xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    xlabel("\$T\$ (\$K\$)")
    xlim(282.72,311.72)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-3)
    ylim(-5,38)
    legend(loc="upper right",numpoints=1,prop={"size"=>legend_font_size})


subplot(223)
    # helper functions/values
    Ts = 286.5:.1:288.5
    T_0 = 287
    T_1 = 287.2
    T_2 = 288
    ΔT_obs = T_2 - T_1

    F_1 = .5
    F_2 = 2
    ΔF = F_2 - F_1

    tw = .02
    a = .058
    N(F,T) = F - 1.28 * (T - T_0) + a * (T - T_0)^2

    ΔQ = N(F_2,T_2) - N(F_1,T_1)
    linear_N(F,T) = F + (-(ΔF - ΔQ) / ΔT_obs) * (T - T_0)


    # legend labels
    plot([0,0],[1,1],label="\$\\lambda = -1.28\$ \$W/m^2/K\$", linestyle="none")
    plot([0,0],[1,1],label="\$a_H = .058\$ \$W/m^2/K^2\$", linestyle="none")
    # curves of N(T)
    plot(Ts,[linear_N(F_1,T) for T=Ts],"k--",color="k")
    plot(Ts,[linear_N(F_2,T) for T=Ts],"k--",linewidth=2,color="k",label="linear")
    plot(Ts,[linear_N(3.71,T) for T=Ts],"k--",linewidth=2,color="0.5")
    plot(Ts,[N(F_1,T) for T=Ts],linewidth=2,color="k")
    plot(Ts,[N(F_2,T) for T=Ts],linewidth=2,color="k",label="quad")
    plot(Ts,[N(3.71,T) for T=Ts],"k",linewidth=2,color="0.5")

    # x axis
    plot(Ts,[0 for T=Ts],"k",linewidth=2)

    # lambda
    plot([T_0-.1,T_0+.1],[N(2,T_0-.1),N(2,T_0+.1)],linewidth=4,color="magenta")
    text(287.02,2.05,"\$\\lambda\$",color="magenta")

    # Ts
    plot([T_0,T_0],[-.08,.08],linewidth=2,color="k")
    plot([T_1,T_1],[-.08,.08],linewidth=2,color="k")
    plot([T_2,T_2],[-.08,.08],linewidth=2,color="k")
    text(T_0-.07,-.4,"\$T_0\$",color="k")
    text(T_1-.07,-.4,"\$T_1\$",color="b")
    text(T_2-.07,-.4,"\$T_2\$",color="r")

    arrow(T_1,N(2,T_2),T_2-T_1-.07,0,shape="full",edgecolor="purple",facecolor="purple",head_width=.08,head_length=.04)
    text(.5 * (T_2+T_1) - .08,N(2,T_2) - .3,"\$\\Delta T\$",color="purple")

    # Qs
    text(T_1-.15,.08,"\$Q_1\$",color="b")
    plot([T_1,T_1],[0,N(F_1,T_1)],"b",linewidth=2)
    plot([T_1-tw,T_1+tw],[0,0],linewidth=2,color="b")
    plot([T_1-tw,T_1+tw],[N(F_1,T_1),N(F_1,T_1)],linewidth=2,color="b")

    text(T_2-.17,.35,"\$Q_2\$",color="r")
    plot([T_2,T_2],[0,N(F_2,T_2)],"r",linewidth=2)
    plot([T_2-tw,T_2+tw],[0,0],linewidth=2,color="r")
    plot([T_2-tw,T_2+tw],[N(F_2,T_2),N(F_2,T_2)],linewidth=2,color="r")

    text(287.22,.42,"\$\\Delta Q\$",color="orange")
    plot([T_1+tw/2,T_1+tw/2],[N(F_1,T_1),N(F_2,T_2)],color="orange",linewidth=2)
    plot([T_1,T_1+tw],[N(F_1,T_1),N(F_1,T_1)],color="orange",linewidth=2)
    plot([T_1,T_1+tw],[N(F_2,T_2),N(F_2,T_2)],color="orange",linewidth=2)

    # Fs
    text(286.85,.18,"\$F_1\$",color="b")
    plot([T_0-tw/2,T_0-tw/2],[0,F_1],"b",linewidth=2)
    plot([T_0-tw,T_0],[0,0],"b",linewidth=2)
    plot([T_0-3tw/2,T_0+tw/2],[F_1,F_1],"b",linewidth=2)

    text(287.015,.9,"\$F_2\$",color="r")
    plot([T_0+tw/2,T_0+tw/2],[0,F_2],"r",linewidth=2)
    plot([T_0,T_0+tw],[0,0],"r",linewidth=2)
    plot([T_0-tw/2,T_0+3tw/2],[F_2,F_2],"r",linewidth=2)

    text(287.2,1.12284,"\$\\Delta F\$",color="g")
    plot([T_1-tw/2,T_1-tw/2],[N(F_1,T_1),N(F_2,T_1)],"g",linewidth=2)
    plot([T_1-tw,T_1],[N(F_1,T_1),N(F_1,T_1)],"g",linewidth=2)
    plot([T_1-3tw/2,T_1+tw/2],[N(F_2,T_1),N(F_2,T_1)],"g",linewidth=2)

    # labels and axes
    annotate("c",xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    xlabel("\$T\$ (\$K\$)")
    xlim(286.5,288.5)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    ylim(-.5,3.5)
    legend(loc="upper right",prop={"size"=>legend_font_size})

subplot(224)
    # helper functions/values
    Ts = 286:.1:300
    T_0 = 287
    T_1 = 287.4
    T_2 = 288
    tw = .02

    # x axis
    plot(Ts,[0 for T=Ts],"k",linewidth=2)

    # legend labels
    plot([0,0],[1,1],label="\$\\lambda = -1.28\$ \$W/m^2/K\$", linestyle="none")
    plot([0,0],[1,1],label="\$a_H = .058\$ \$W/m^2/K^2\$", linestyle="none")
    # curves of N(T)
    plot(Ts,[linear_N(F_1,T) for T=Ts],"k--",linewidth=2,color="0.8")
    plot(Ts,[N(F_1,T) for T=Ts],"k",linewidth=2,color="0.8")
    plot(Ts,[linear_N(2,T) for T=Ts],"k--",linewidth=2,color="0.8")
    plot(Ts,[N(2,T) for T=Ts],"k",linewidth=2,color="0.8")
    plot(Ts,[linear_N(3.71,T) for T=Ts],"b--",linewidth=2,label="\$2xCO_2\$ (linear), \$\\Delta T_{2x} = $(round(ΔT(3.71,-1.28,0),1))K\$")
    plot(Ts,[N(3.71,T) for T=Ts],"b",linewidth=2,label="\$2xCO_2\$ (quad), \$\\Delta T_{2x} = $(round(ΔT(-1.28,.058,3.71),1))K\$")
    plot(Ts,[linear_N(7.42,T) for T=Ts],"r--",linewidth=2,label="\$4xCO_2\$ (linear), \$\\Delta T_{4x} = $(round(ΔT(7.42,-1.28,0),1))K\$")
    plot(Ts,[N(7.42,T) for T=Ts],"r",linewidth=2,label="\$4xCO_2\$ (quad), \$\\Delta T_{4x} = \$ ?")

    # zoom box
    plot([286.5,288.5],[-1.,-1.],"k",linewidth=2)
    plot([286.5,286.5],[-1.,3.],"k",linewidth=2)
    plot([286.5,288.5],[3.,3.],"k",linewidth=2)
    plot([288.5,288.5],[-1.,3.],"k",linewidth=2)
    plot(286.5:.01:288.24,[N(F_1,T) for T=286.5:.01:288.24],"k",linewidth=2,color="0.5")
    plot(286.5:.1:288.5,[N(2,T) for T=286.5:.1:288.5],"k",linewidth=2,color="0.5")

    # label
    plot([T_0-.1,T_0+.1],[N(2,T_0-.1),N(2,T_0+.1)],linewidth=4,color="magenta")

    # tick
    text(T_0-.4,-2.5,"\$T_0\$",color="k")
    plot([T_0,T_0],[-.4,.4],linewidth=2,color="k")

    # labels and axes
    annotate("d",xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    xlabel("\$T\$ (\$K\$)")
    xlim(286,300)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    ylim(-7,21)
    legend(loc="upper right",prop={"size"=>legend_font_size})

subplots_adjust(left=.1,right=.97,wspace=.2,hspace=.25,bottom=.1,top=.95)
savefig("figs/1.eps")

# ----------------
# Figure 2
# ----------------
figure(2,(15,5))

# GCM data
models = [
    {
        :name => "Stouffer and Manabe\n(2003) - GFDL",
        :T_s => [-7.85, 0.0, 4.33, 7.9],
        :Net => [3.71,0.0,-3.71,-7.42],
        :T_2x => 4.33,
        :T_4x => 7.9,
        :T_8x => 10000 # value unknown
    },
    {
        :name => "Hansen et al.\n(2005) - GISS",
        :T_s => [-4.56, -3.32, -1.78, .0, 0.58, 1.10, 1.96, 4.06, 7.02],
        :Net => [12.68,8.88,4.61,-0.0,-1.44,-2.64,-4.52,-9.27,-14.65],
        :T_2x => 1.96,
        :T_4x => 4.06,
        :T_8x => 7.02
    },
    {
        :name => "Colman and McAvaney\n(2009) - BMRC",
        :T_s => [-12.3875, -7.74805, -4.60923, -2.70842, 0., 2.49182, 5.11104, 6.84392, 8.9713, 12.7421],
        :Net => [15.392,11.704,7.936,4.04,-0.0,-4.28,-8.8,-13.776,-19.848,-26.24],
        :T_2x => 2.49182,
        :T_4x => 5.11104,
        :T_8x => 6.84392
    },
    {
        :name => "Jonko et al.\n(2013) - CAM-SOM",
        :T_s => [0,2.2,4.7,7.8],
        :Net => [0.0,-3.61,-7.52,-11.29],
        :T_2x => 2.2,
        :T_4x => 4.7,
        :T_8x => 7.8
    },
    {
        :name => "Meraner et al.\n(2013) - ECHAM6",
        :T_s => [0,2.79,6.7,12.46,22.68],
        :Net => [0,-4.26,-9.09,-14.24,-19.18],
        :T_2x => 2.79,
        :T_4x => 6.7,
        :T_8x => 12.46
    }
]

# helper functions/values
function polyfit_exps(x,y,n,flag=false)
    A = float([ float(x[i])^p for i = 1:length(x), p = (flag == :no_constant ? 1 : 0):n ])
    A \ y
end

function run_model(m,F)
    markersize = (m[:name][1:2] == "Me" ? 21 : (m[:name][1] == 'H' ? 17 : 15))
    N0, λ, a = [round(x,3) for x=polyfit_exps(m[:T_s], m[:Net],2)]
    # λ, a = [round(x,3) for x=polyfit_exps(m[:T_s], m[:Net],2,:no_constant)]
    plot(a,(F==3.71? m[:T_2x] : (F==7.42 ? m[:T_4x] : m[:T_8x]) ),marker="\$$(m[:name][1])\$",markersize=markersize,linestyle="none",color="k",label=m[:name])
end

# range of λs from CMIP5, assuming linearity (Andrews et al., 2012)
CMIP5_range = {
    :max => [-.79, "purple"],
    :mean => [-1.17, "green"],
    :min => [-1.78, "orange"]
}

as = -.1:.000005:.1
pos_as = .001:.001:.1
neg_as = -.1:.001:-.001

function subplot_setup(i)
    # helper values and functions
    doublings = 2^i
    F = 3.71i
    y_max = 10i
    function plot_pair(kv_pair,F)
        stat, pair = kv_pair
        λ, color = pair
        plot(0,-F/λ,marker=".",color=color,markersize=20,linestyle="none")
        plot(as,[ΔT(F,λ,a) for a=as], color=color, label="\$\\lambda_{$stat} \$ = $λ", linewidth=2)
    end

    # plot GCMs
    plot(1000,1000,linestyle="none",label="GCM experiments:")
    [run_model(m,F) for m=models]
    # plot lines
    plot(1000,1000,linestyle="none",label=" ")
    plot(1000,1000,linestyle="none",label="CMIP5 \$\\lambda\$s:")
    [plot_pair(pair,F) for pair=CMIP5_range]

    # unstable region
    fill_between(pos_as,[sqrt(F/a) for a=pos_as],30,color="0.9",facecolor="0.9")

    # disallowed region
    fill_between(neg_as,[sqrt(-F/a) for a=neg_as],30,color="0.7",facecolor="0.7")
    plot(1000,1000,linestyle="none",label=" ")
    plot(neg_as,[sqrt(-F/a) for a=neg_as], color="black", linewidth=2,label="\$\\lambda\$ = 0 \$W/m^2/K\$")

    # center line
    plot([0 for f=[0,y_max]],[0,y_max],color="0.4")

    # labels and axes
    annotate(string("abc"[i]),xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    title("\$$(doublings)xCO_2\$")
    xlabel("\$a\$ (\$W/m^2/K^2\$)")
    xlim(-.1,.1)
    ylabel("\$\\Delta T_{$(doublings)x}\$ (\$K\$)",labelpad=-5)
    ylim(0,y_max)
end

subplot(131)
    text(.0535,8.7,"unstable",fontsize=13)

    # Figure 1c,d case
    plot(.058,ΔT(3.71,-1.28,.058),"k",marker="*",markersize=16)

    subplot_setup(1)

subplot(132)
    text(.049,14.8,"unstable",fontsize=13)

    # CMIP5 constraint case
    spec_as = .001:.0001:.1
    # RCP8.5 goes up to 1962ppm
    F_limit = log2(1962/270) * 3.71
    derived_λ(F,a,ΔT) = -(F+a*ΔT^2)/ΔT
    CMIP5_limit(a) = ΔT(2*3.71,derived_λ(F_limit,a,sqrt(F_limit/a)),a)
    # draw region between min ΔT_{4x}, CMIP5 limit
    fill_between(spec_as,[-2*3.71/CMIP5_range[:min][1] for a=spec_as],[min(CMIP5_limit(a),F_limit) for a=spec_as],facecolor="1.0",color="pink")

    subplot_setup(2)


subplot(133)
    text(.049,14.8*1.5,"unstable",fontsize=13)

    subplot_setup(3)
    legend(loc="center left", bbox_to_anchor=(1.1,.43),numpoints=1,prop={"size"=>13},frameon=false)


subplots_adjust(left=.06,right=.8,bottom=.2,top=.85,wspace=.26)
savefig("figs/2.eps")

# ----------------
# Figure 3
# ----------------
figure(3,(15,8))

# y pos shows order in which display in bottom subplot
ΔT_4x_dists = {
    :linear => {:color => "black", :linestyle => "--", :y_pos => 2, :label => "linear"},
    :a_min => {:color => "blue", :linestyle => "-", :y_pos => 1, :label => "\$a_C=-.035\$"},
    :a_max => {:color => "red", :linestyle => "-", :y_pos => 5, :label => "\$a_H=.058\$\nP(unstable)=.17,.55"},
    :a_dist => {:color => "purple", :linestyle => "-", :y_pos => 3, :label => "\$a \\sim unif(\\pm .06)\$\nP(unstable)=.03,.14"},
    :a_pos_dist => {:color => "green", :linestyle => "-", :y_pos => 4, :label => "\$a \\sim unif(0, .06)\$\nP(unstable)=.07,.28"}
}

function column_setup(i)
    epdf(cdf,step_size=.1) = ΔT -> (-cdf(ΔT+2step_size) + 8*cdf(ΔT+step_size) - 8*cdf(ΔT-step_size)+cdf(ΔT-2step_size)) / (12step_size)

    for (k,v) in ΔT_4x_dists
        v[:cdf] = ecdf(v[:values])
        v[:pdf] = epdf(v[:cdf])
    end

    ΔT_range = 0:.1:15

    # PDFs
    subplot2grid((11,2),(0,i-1),rowspan=4)
        title(i == 1 ? "Otto et al. (2013)" : "Murphy et al. (2009)", fontsize=16)

        # curves
        for (k,v) in ΔT_4x_dists
            plot(ΔT_range, [v[:pdf](ΔT) for ΔT=ΔT_range], color=v[:color], linestyle=v[:linestyle], linewidth=2)
        end

        # labels and axes
        annotate(string("ad"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.88))
        xlim(0,15)
        ylabel("probability density (K\$^{-1}\$)", fontsize=16)
        ylim(0,.4)
        yticks(0:.1:.4,0:.1:.4)

    # CDFs
    subplot2grid((33,2),(14,i-1),rowspan=12)

        # curves
        for (k,v) in ΔT_4x_dists
            plot(ΔT_range, [v[:cdf](ΔT) for ΔT=ΔT_range], color=v[:color], linestyle=v[:linestyle], linewidth=2)
        end

        # labels and axes
        annotate(string("be"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.88))
        xlim(0,15)
        ylabel("cumulative probability", fontsize=16)

    # percentiles
    subplot2grid((6,2),(5,i-1))

        ΔT_4x_limits(ΔT_4x_cdf) = [try fzeros(ΔT_4x -> ΔT_4x_cdf(ΔT_4x) - perc, 0, 20)[1] catch 1e6 end for perc=[.05,.17,.5,.83,.95]]

        msize=9
        for k in [:linear,:a_min,:a_max,:a_dist,:a_pos_dist]
            v = ΔT_4x_dists[k]
            plot(ΔT_4x_limits(v[:cdf]),[.0125 - .025v[:y_pos] for i=1:5],color=v[:color],linestyle=v[:linestyle],marker="|",markersize=msize)
            plot(ΔT_4x_limits(v[:cdf])[2:4],[.0125 - .025v[:y_pos] for i=1:3],color=v[:color],linewidth=3)
            plot([1000,1000],[1000,1000],color=v[:color],linestyle=v[:linestyle],linewidth=3,label=v[:label])
        end
        cmip_ΔTs = 2*[3.69, 3.25, 4.08, 3.97, 2.39, 2.44, 4.59, 2.08, 4.13, 4.67, 2.72, 3.63, 3.45, 2.60, 2.80]
        plot(cmip_ΔTs,[.01 for c=cmip_ΔTs],"kx",markersize=10,label="CMIP5 models")

        annotate(string("cf"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.75))
        xlabel("\$\\Delta T_{4x}\$ (K)")
        xlim(0,15)
        ylim(-.14,.03)
        yticks([],[])
end

samples = int(5e6)


    # Otto et al. (2013)
    ΔQs = rand(Normal(0.65,.164148),samples)
    ΔTs = rand(Normal(0.75,.12159),samples)
    ΔF_ghgs = rand(Normal(2.83,.17),samples)
    ΔF_not_ghgs = rand(Normal(-.88,.31),samples)
    ΔFs = ΔF_ghgs + ΔF_not_ghgs
    ΔF_4xs = ΔF_ghgs * (3.44/2.83) * 2

    a_min = -.035
    a_max = .058
    a_dist = rand(Uniform(-.06,.06),samples)
    a_pos_dist = rand(Uniform(0.,.06),samples)

    λ_linears = (ΔQs - ΔFs) ./ ΔTs
    λ_a_mins = (ΔQs - ΔFs - a_min * ΔTs.^2) ./ ΔTs
    λ_a_maxs = (ΔQs - ΔFs - a_max * ΔTs.^2) ./ ΔTs
    λ_a_dists = (ΔQs - ΔFs - a_dist .* ΔTs.^2) ./ ΔTs
    λ_a_pos_dists = (ΔQs - ΔFs - a_pos_dist .* ΔTs.^2) ./ ΔTs

    ΔT_4x_dists[:linear][:values] = ΔF_4xs ./ -λ_linears
    ΔT_4x_dists[:a_min][:values] = Float64[ΔT(ΔF_4xs[i],λ_a_mins[i],a_min) for i=1:samples]
    ΔT_4x_dists[:a_max][:values] = Float64[ΔT(ΔF_4xs[i],λ_a_maxs[i],a_max) for i=1:samples]
    ΔT_4x_dists[:a_dist][:values] = Float64[ΔT(ΔF_4xs[i],λ_a_dists[i],a_dist[i]) for i=1:samples]
    ΔT_4x_dists[:a_pos_dist][:values] = Float64[ΔT(ΔF_4xs[i],λ_a_pos_dists[i],a_pos_dist[i]) for i=1:samples]

    column_setup(1)

    # Murphy et al. (2009)
    ΔF_4x = 7.42

    λ_dist = Normal(-1.25,.5)
    λs = rand(λ_dist,samples)

    ΔT_4x_dists[:linear][:values] = ΔF_4x ./ -λs
    ΔT_4x_dists[:a_min][:values] = Float64[ΔT(ΔF_4x,λs[i],a_min) for i=1:samples]
    ΔT_4x_dists[:a_max][:values] = Float64[ΔT(ΔF_4x,λs[i],a_max) for i=1:samples]
    ΔT_4x_dists[:a_dist][:values] = Float64[ΔT(ΔF_4x,λs[i],a_dist[i]) for i=1:samples]
    ΔT_4x_dists[:a_pos_dist][:values] = Float64[ΔT(ΔF_4x,λs[i],a_pos_dist[i]) for i=1:samples]

    column_setup(2)
    legend(loc="center left",numpoints=1,prop={"size"=>16},frameon=false,bbox_to_anchor=(1,3.5))

subplots_adjust(right=.8)
savefig("figs/3.eps")


# ----------------
# Lewis and Curry recalculation
# ----------------
# 1) unzip GMST.zip
# 2) run lewis_and_curry.r to generate csvs
# 3) uncomment what follows:
# println("Lewis and Curry recalculations:")
# a_min_cdf = ecdf(readtable("data/lewis_and_curry_a_min.csv")[:x_1])
# println("P(ΔT_4x > 8K), a = a_min: $(int(round(1-a_min_cdf(8),2)*100))%")
# a_max_cdf = ecdf(readtable("data/lewis_and_curry_a_max.csv")[:x_1])
# println("P(ΔT_4x > 8K), a = a_max: $(int(round(1-a_max_cdf(8),2)*100))%")
# a_pos_dist_cdf = ecdf(readtable("data/lewis_and_curry_a_pos_dist.csv")[:x_1])
# println("P(unstable), a ~ unif(0,.06): $(int(round(1-a_pos_dist_cdf(1e3),2)*100))%")
