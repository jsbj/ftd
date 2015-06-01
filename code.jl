# code for "Feedback Temperature Dependence and Equilibrium Climate Sensitivity"
# Jonah Bloch-Johnson, Raymond Pierrehumbert, Dorian Abbot
# this code is written in Julia (http://julialang.org/)
using PyPlot, Roots, DataFrames, Distributions

plot_these = 4:5 # 1:6

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

# warming for quintic model
function ΔT(F,λ,a,g)
    try
        sort(fzeros(ΔT -> F + λ * ΔT + a * ΔT^2 + g * ΔT^5, 0, 50))[1]
    catch
        Inf
    end
end

# warming for quintic model
function ΔT_w_CO2(D,λ,a,b)
    F = 3.71D
    λ = λ + D * b
    if a == 0
        -F / λ
    else
        dtrx = λ^2 - 4a * F
        dtrx >= 0. ? (-λ - sqrt(dtrx))/2a : Inf
    end
end


subplot_label_x = .025
subplot_label_y = .91
legend_font_size = 12

# ----------------
# Figure 1
# ----------------
if in(1, plot_these)
    figure(1,(10,8))

    subplot(221)
        # helper functions/values
        λ = -.88
        N(T,a) = 3.71 + λ * (T - 287) + a * (T - 287)^2
        Ts = 285:.2:297

        # legend label
        plot([0,0],[1,1],label="\$\\lambda = $λ\$ \$W/m^2/K\$",linestyle="none")
        # curves of N(T)
        plot(Ts,[N(T,    0) for T=Ts],"k--",label="linear, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,0),1))K\$",linewidth=2)
        plot(Ts,[N(T,-.035) for T=Ts],"b",linewidth=2,label="\$a_C=-0.035\$, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,-.035),1))K\$")
        plot(Ts,[N(T, .03) for T=Ts],"g",linewidth=2,label="\$a_M=0.03\$, \$\\Delta T_{2x} = $(round(ΔT(3.71,λ,.03),1))K\$")
        plot(Ts,[N(T, .058) for T=Ts],"r",linewidth=2,label="\$a_H=0.058\$, \$\\Delta T_{2x} = \$ ?")

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
        plot([1000,2000],[1000,1000],linestyle="none",label="ECHAM6, 32xCO\$_2\$")
        # curves/plots of N(T)
        plot(Ts,Ns,"bo",label="\$\\Delta T\$ vs. \$F_{32x} - F\$")
        plot(T_grid+287.72,[27 - 1.52 * T for T=T_grid],"k--",linewidth=2,label="linear (\$\\lambda_M = -1.52\$)")
        plot(T_grid+287.72,[27 - 1.52 * T + .03 * T^2 for T=T_grid],"r-",linewidth=2,label="quad (\$a_M = 0.03\$)")
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
        plot([0,0],[1,1],label="\$a_H = 0.058\$ \$W/m^2/K^2\$", linestyle="none")
        # curves of N(T)
        plot(Ts,[linear_N(F_1,T) for T=Ts],"k--",color="k")
        plot(Ts,[linear_N(F_2,T) for T=Ts],"k--",linewidth=2,color="k",label="linear")
        plot(Ts,[linear_N(3.71,T) for T=Ts],"k--",linewidth=2,color="blue")
        plot(Ts,[N(F_1,T) for T=Ts],linewidth=2,color="k")
        plot(Ts,[N(F_2,T) for T=Ts],linewidth=2,color="k",label="quad")
        plot(Ts,[N(3.71,T) for T=Ts],"k",linewidth=2,color="blue")

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
        plot([0,0],[1,1],label="\$a_H = 0.058\$ \$W/m^2/K^2\$", linestyle="none")
        # curves of N(T)
        plot(Ts,[linear_N(F_1,T) for T=Ts],"k--",linewidth=2,color="k")
        plot(Ts,[N(F_1,T) for T=Ts],"k",linewidth=2,color="k")
        plot(Ts,[linear_N(2,T) for T=Ts],"k--",linewidth=2,color="k")
        plot(Ts,[N(2,T) for T=Ts],"k",linewidth=2,color="k")
        plot(Ts,[linear_N(3.71,T) for T=Ts],"b--",linewidth=2,label="\$2xCO_2\$ (linear), \$\\Delta T_{2x} = $(round(ΔT(3.71,-1.28,0),1))K\$")
        plot(Ts,[N(3.71,T) for T=Ts],"b",linewidth=2,label="\$2xCO_2\$ (quad), \$\\Delta T_{2x} = $(round(ΔT(3.71,-1.28,.058),1))K\$")
        plot(Ts,[linear_N(7.42,T) for T=Ts],"r--",linewidth=2,label="\$4xCO_2\$ (linear), \$\\Delta T_{4x} = $(round(ΔT(7.42,-1.28,0),1))K\$")
        plot(Ts,[N(7.42,T) for T=Ts],"r",linewidth=2,label="\$4xCO_2\$ (quad), \$\\Delta T_{4x} = \$ ?")

        # zoom box
        plot([286.5,288.5],[-1.,-1.],"k",linewidth=2)
        plot([286.5,286.5],[-1.,3.],"k",linewidth=2)
        plot([286.5,288.5],[3.,3.],"k",linewidth=2)
        plot([288.5,288.5],[-1.,3.],"k",linewidth=2)
        # plot(286.5:.01:288.24,[N(F_1,T) for T=286.5:.01:288.24],"k",linewidth=2,color="0.5")
        # plot(286.5:.1:288.5,[N(2,T) for T=286.5:.1:288.5],"k",linewidth=2,color="0.5")

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
end

# ----------------
# Figure 2
# ----------------
if in(2,plot_these)
    figure(2,(15,5))

    # GCM data
    models = [
        {
            :name => "Stouffer and Manabe\n(2003) - GFDL",
            :T_s => [-7.85, 0.0, 4.33, 7.9],
            :Net => [3.71,0.0,-3.71,-7.42],
            :D => [-1,0,1,2],
            :T_2x => 4.33,
            :T_4x => 7.9,
            :T_8x => 10000 # value unknown
        },
        {
            :name => "Hansen et al.\n(2005) - GISS-E",
            :T_s => [-4.56, -3.32, -1.78, .0, 0.58, 1.10, 1.96, 4.06, 7.02],
            :Net => [12.68,8.88,4.61,-0.0,-1.44,-2.64,-4.52,-9.27,-14.65],
            :D => [-3,-2,-1,0,.322,.585,1,2,3],
            :T_2x => 1.96,
            :T_4x => 4.06,
            :T_8x => 7.02
        },
        {
            :name => "Colman and McAvaney\n(2009) - BMRC",
            :T_s => [-12.3875, -7.74805, -4.60923, -2.70842, 0., 2.49182, 5.11104, 6.84392, 8.9713, 12.7421],
            :Net => [15.392,11.704,7.936,4.04,-0.0,-4.28,-8.8,-13.776,-19.848,-26.24],
            :D => [-4,-3,-2,-1,0,1,2,3,4,5],
            :T_2x => 2.49182,
            :T_4x => 5.11104,
            :T_8x => 6.84392
        },
        {
            :name => "Jonko et al.\n(2013) - CAM3",
            :T_s => [0,2.2,4.7,7.8],
            :Net => [0.0,-3.61,-7.52,-11.29],
            :D => [0,1,2,3],
            :T_2x => 2.2,
            :T_4x => 4.7,
            :T_8x => 7.8
        },
        {
            :name => "Meraner et al.\n(2013) - ECHAM6",
            :T_s => [0,2.79,6.7,12.46,22.68],
            :Net => [0,-4.26,-9.09,-14.24,-19.18],
            :D => [0,1,2,3,4],
            :T_2x => 2.79,
            :T_4x => 6.7,
            :T_8x => 12.46
        }
    ]
    #
    # # helper functions/values
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
        plot(1000,1000,linestyle="none",label="GCMs:")
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
        phrase = "runaway"
        text(.0535,8.7,phrase,fontsize=13)

        # Figure 1c,d case
        # plot(.058,ΔT(3.71,-1.28,.058),"k",marker="*",markersize=16)

        subplot_setup(1)

    subplot(132)
        text(.049,14.8,phrase,fontsize=13)

        # CMIP5 constraint case
        spec_as = -.1:.0001:.1
        # RCP8.5 goes up to 1962ppm
        F_limit = log2(1962/270) * 3.71
        derived_λ(F,a) = -(2*sqrt(F*a))
        # derived_λ(F,a,ΔT) = -(F+a*ΔT^2)/ΔT
        CMIP5_limit(a) = a > 0 ? ΔT(2*3.71,derived_λ(F_limit,a),a) : (a < 0 ? sqrt(2*3.71 / -a) : Inf)
        # CMIP5_limit(a) = a > 0 ? ΔT(2*3.71,derived_λ(F_limit,a,sqrt(F_limit/a)),a) : (a < 0 ? sqrt(2*3.71 / -a) : Inf)
        # draw region between min ΔT_{4x}, CMIP5 limit
        fill_between(spec_as,[-2*3.71/CMIP5_range[:min][1] for a=spec_as],[min(CMIP5_limit(a),F_limit) for a=spec_as],facecolor="#FFF0F0",color="#FFF0F0")

        subplot_setup(2)


    subplot(133)
        text(.049,14.8*1.5,phrase,fontsize=13)

        subplot_setup(3)
        legend(loc="center left", bbox_to_anchor=(1.1,.43),numpoints=1,prop={"size"=>13},frameon=false)


    subplots_adjust(left=.06,right=.8,bottom=.2,top=.85,wspace=.26)
    savefig("figs/2.eps")
end

# ----------------
# Figure 3
# ----------------
if in(3,plot_these)
    figure(3,(15,8))

    # y pos shows order in which display in bottom subplot
    ΔT_4x_dists = {
        :linear => {:color => "black", :linestyle => "--", :y_pos => 2, :label => "linear"},
        :a_min => {:color => "blue", :linestyle => "-", :y_pos => 1, :label => "\$a_C=-0.035\$"},
        :a_max => {:color => "red", :linestyle => "-", :y_pos => 5, :label => "\$a_H=0.058\$\n\$P_{qr}=11\\%,17\\%,55\\%\$"},
        :a_dist => {:color => "purple", :linestyle => "-", :y_pos => 3, :label => "\$a \\sim unif(\\pm 0.06)\$\n\$P_{qr}=3\\%,3\\%,14\\%\$"},
        :a_pos_dist => {:color => "green", :linestyle => "-", :y_pos => 4, :label => "\$a \\sim unif(0, 0.06)\$\n\$P_{qr}=6\\%,7\\%,28\\%\$"}
    }
    
    
    function calculate_curves(i, study, ΔT_range = 0:.1:15)
        for case in keys(ΔT_4x_dists)
            ΔT_4x_dists[case][:values] = readtable("data/$(study)_$case.csv")[:x_1]  
        end
        
        epdf(cdf,step_size=.1) = ΔT -> (-cdf(ΔT+2step_size) + 8*cdf(ΔT+step_size) - 8*cdf(ΔT-step_size)+cdf(ΔT-2step_size)) / (12step_size)

        for (k,v) in ΔT_4x_dists
            v[:cdf] = ecdf(v[:values])
            v[:pdf] = epdf(v[:cdf])
        end
    
        curves = DataFrame()
        
        for (k,v) in ΔT_4x_dists
            curves[symbol("$(k)_pdf")] = [round(v[:pdf](ΔT),6) for ΔT=ΔT_range]
            curves[symbol("$(k)_cdf")] = [round(v[:cdf](ΔT),6) for ΔT=ΔT_range]
        end
    
        writetable("data/$(study)_curves.csv", curves)
    end

    function calculate_percentiles(i, study)
        percentiles = DataFrame()
        
        ΔT_percentiles(ΔT_4x_cdf) = [try fzeros(ΔT_4x -> ΔT_4x_cdf(ΔT_4x) - perc, 0, 20)[1] catch 1e6 end for perc=[.05,.17,.5,.83,.95]]
        
        for (k,v) in ΔT_4x_dists
            percentiles[k] = ΔT_percentiles(ecdf(readtable("data/$(study)_$k.csv")[:x_1]))
        end
    
        writetable("data/$(study)_percentiles.csv", percentiles)
    end

    function column_setup(i, study)
        
        curves = readtable("data/$(study)_curves.csv")

        ΔT_range = 0:.1:15

        # PDFs
        subplot2grid((11,3),(0,i-1),rowspan=4)
            title(["Lewis and Curry (2014)", "Otto et al. (2013)", "Murphy et al. (2009)"][i], fontsize=16, fontweight="light")

            # curves
            for (k,v) in ΔT_4x_dists
                plot(ΔT_range, curves[symbol("$(k)_pdf")], color=v[:color], linestyle=v[:linestyle], linewidth=2)
            end

            # labels and axes
            annotate(string("adg"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.88))
            xlim(0,15)
            if i == 1
                ylabel("probability density (K\$^{-1}\$)", fontsize=16, fontweight="light")
            end
            ylim(0,.5)
            yticks(0:.1:.5,0:.1:.5)

        # CDFs
        subplot2grid((33,3),(14,i-1),rowspan=12)

            # curves
            for (k,v) in ΔT_4x_dists
                plot(ΔT_range, curves[symbol("$(k)_cdf")], color=v[:color], linestyle=v[:linestyle], linewidth=2)
            end

            # labels and axes
            annotate(string("beh"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.88))
            xlim(0,15)
            if i == 1
                ylabel("cumulative probability", fontsize=16, fontweight="light")
            end

        # percentiles
        subplot2grid((6,3),(5,i-1))
        
            percentiles = readtable("data/$(study)_percentiles.csv")

            msize=9
            for k in [:a_min,:linear,:a_dist,:a_pos_dist,:a_max]
                v = ΔT_4x_dists[k]
                plot(percentiles[k],[.0125 - .025v[:y_pos] for i=1:5],color=v[:color],linestyle=v[:linestyle],marker="|",markersize=msize)
                plot(percentiles[k][2:4],[.0125 - .025v[:y_pos] for i=1:3],color=v[:color],linewidth=3)
                plot([1000,1000],[1000,1000],color=v[:color],linestyle=v[:linestyle],linewidth=3,label=v[:label])
            end
            cmip_ΔTs = 2*[3.69, 3.25, 4.08, 3.97, 2.39, 2.44, 4.59, 2.08, 4.13, 4.67, 2.72, 3.63, 3.45, 2.60, 2.80]
            plot(cmip_ΔTs,[.01 for c=cmip_ΔTs],"kx",markersize=10,label="CMIP5 models\n(abrupt4xCO\$_2\$ runs)")

            annotate(string("cfi"[i]),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom",xy=(.02,.75))
            xlabel("\$\\Delta T_{4x}\$ (K)")
            xlim(0,15)
            ylim(-.14,.03)
            yticks([],[])
    end

    samples = int(5e6)

    # Lewis and Curry (2014)
    # ----------------------

    # to calculate
    # 1) unzip GMST.zip
    # 2) run lewis_and_curry.r to generate csvs
    # calculate_curves(1, "lewis_and_curry")
    # calculate_percentiles(1, "lewis_and_curry")

    column_setup(1, "lewis_and_curry")

    #
    # # Otto et al. (2013)
    # # ------------------
    #
    # # to calculate
    # ΔQs = rand(Normal(0.65,.164148),samples)
    # ΔTs = rand(Normal(0.75,.12159),samples)
    # ΔF_ghgs = rand(Normal(2.83,.17),samples)
    # ΔF_not_ghgs = rand(Normal(-.88,.31),samples)
    # ΔFs = ΔF_ghgs + ΔF_not_ghgs
    # ΔF_4xs = ΔF_ghgs * (3.44/2.83) * 2
    #
    # a_min = -.035
    # a_max = .058
    # a_dist = rand(Uniform(-.06,.06),samples)
    # a_pos_dist = rand(Uniform(0.,.06),samples)
    #
    # λ_linears = (ΔQs - ΔFs) ./ ΔTs
    # λ_a_mins = (ΔQs - ΔFs - a_min * ΔTs.^2) ./ ΔTs
    # λ_a_maxs = (ΔQs - ΔFs - a_max * ΔTs.^2) ./ ΔTs
    # λ_a_dists = (ΔQs - ΔFs - a_dist .* ΔTs.^2) ./ ΔTs
    # λ_a_pos_dists = (ΔQs - ΔFs - a_pos_dist .* ΔTs.^2) ./ ΔTs
    #
    #
    # study = "otto_et_al"
    # writetable("data/$(study)_linear.csv", DataFrame(x_1 = ΔF_4xs ./ -λ_linears))
    # writetable("data/$(study)_a_min.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4xs[i],λ_a_mins[i],a_min) for i=1:samples]))
    # writetable("data/$(study)_a_max.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4xs[i],λ_a_maxs[i],a_max) for i=1:samples]))
    # writetable("data/$(study)_a_dist.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4xs[i],λ_a_dists[i],a_dist[i]) for i=1:samples]))
    # writetable("data/$(study)_a_pos_dist.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4xs[i],λ_a_pos_dists[i],a_pos_dist[i]) for i=1:samples]))
    # calculate_curves(2,"otto_et_al")
    # calculate_percentiles(2,"otto_et_al")
    column_setup(2,"otto_et_al")
    
    # # Murphy et al. (2009)
    # # --------------------
    #
    # # to calculate:
    # ΔF_4x = 7.42
    #
    # λ_dist = Normal(-1.25,.5)
    # λs = rand(λ_dist,samples)
    #
    # study = "murphy_et_al"
    # writetable("data/$(study)_linear.csv", DataFrame(x_1 = ΔF_4x ./ -λs))
    # writetable("data/$(study)_a_min.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4x,λs[i],a_min) for i=1:samples]))
    # writetable("data/$(study)_a_max.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4x,λs[i],a_max) for i=1:samples]))
    # writetable("data/$(study)_a_dist.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4x,λs[i],a_dist[i]) for i=1:samples]))
    # writetable("data/$(study)_a_pos_dist.csv", DataFrame(x_1 = Float64[ΔT(ΔF_4x,λs[i],a_pos_dist[i]) for i=1:samples]))
    # calculate_curves(3, "murphy_et_al")
    # calculate_percentiles(3, "murphy_et_al")
    column_setup(3, "murphy_et_al")
    legend(loc="center left",numpoints=1,prop={"size"=>16},frameon=false,bbox_to_anchor=(1,3.5),labelspacing=1)

    subplots_adjust(left=.07,right=.79)
    savefig("figs/3.eps")
end

# ----------------
# Figure 4
# ----------------
if in(4,plot_these)
    figure(4,(15,6))
    Ts = 282:.2:302
    
    λ = -.88
    N(T,F,b=-.0001,g=0,a=.058) = F + λ * (T - 287) + a * (T - 287)^2 + b * (T - 287)^3 + g * (T - 287)^5
    

    subplot2grid((1,4),(0,0))
    a_1 = .04
    annotate(string("a"),xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    plot(Ts, [N(T,4,0,0,0) for T=Ts], linestyle = "--", color="0.5", linewidth=2)
    plot(Ts, [N(T,4,0,0,a_1) for T=Ts], color="0.5", linewidth=2)
    plot(Ts, [N(T,4,0,-4e-6,a_1) for T=Ts], linewidth=2)
    plot(Ts, [N(T,4,-1.3e-3,0,a_1) for T=Ts], linewidth=2)
    plot(Ts, [N(T,4,-7.5e-4,0,a_1) for T=Ts], linewidth=2)
    # plot(Ts, [N(T,F,0,-3e-6) for T=Ts], linewidth=2)
    ylim(-2,15)
    xlim(Ts[1],Ts[end])
    # x axis
    plot(Ts,[0 for T=Ts],"k",linewidth=2)
    plot([287,287],[-.2,.2],"k",linewidth=2)
    xlabel("\$T\$ (\$K\$)",labelpad=-2.5)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    text(286,-1,"\$T_0\$")
    xticks(285:5:300,285:5:300)
    arrow(287,.4,0,3,shape="full",facecolor="k",head_width=.4,head_length=.2)
    arrow(287.7,-.8,4.8,0,shape="full",edgecolor="r",facecolor="r",head_width=.2,head_length=.4)
    arrow(287.7,-.6,4.8,0,shape="full",edgecolor="b",facecolor="b",head_width=.2,head_length=.4)
    arrow(287.7,-.4,4.8,0,shape="full",edgecolor="g",facecolor="g",head_width=.2,head_length=.4)
    arrow(287.7,-.2,3.2,0,shape="full",edgecolor="k",facecolor="k",head_width=.2,head_length=.4)
    text(285.3,1.7,"\$F\$")
    text(287 + 1.2,-1.6,"\$\\Delta T\$")



    subplot2grid((1,4),(0,1),colspan=3)
    annotate(string("b"),xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
    Ts = 282:.2:350
    plot(Ts, [N(T,4,0,0,0) for T=Ts], linestyle = "--", color="0.5", linewidth=2)
    plot(Ts, [N(T,4,0,0) for T=Ts], color="0.5", linewidth=2)
    plot(Ts, [N(T,4,0,-4e-6) for T=Ts], linewidth=2)
    plot(Ts, [N(T,4,-1.3e-3,0) for T=Ts], linewidth=2)
    plot(Ts, [N(T,4,-7.5e-4,0) for T=Ts], linewidth=2)
    # plot(Ts, [N(T,F,0,-3e-6) for T=Ts], linewidth=2)
    ylim(-2,15)
    xlim(Ts[1],Ts[end])
    # x axis
    plot(Ts,[0 for T=Ts],"k",linewidth=2)
    plot([287,287],[-.2,.2],"k",linewidth=2)
    xlabel("\$T\$ (\$K\$)",labelpad=-2.5)
    ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    text(286,-1,"\$T_0\$")

    arrow(287,.4,0,3,shape="full",facecolor="k",head_width=.4,head_length=.2)
    arrow(287.7,-.8,57.7,0,shape="full",edgecolor="r",facecolor="r",head_width=.2,head_length=.4)
    arrow(287.7,-.6,15.7,0,shape="full",edgecolor="b",facecolor="b",head_width=.2,head_length=.4)
    arrow(287.7,-.4,6.2,0,shape="full",edgecolor="g",facecolor="g",head_width=.2,head_length=.4)
    arrow(287.7,-.2,3.2,0,shape="full",edgecolor="k",facecolor="k",head_width=.2,head_length=.4)
    text(285.3,1.7,"\$F\$")
    text(287 + 1.2,-1.6,"\$\\Delta T\$")

    subplots_adjust(left=.06,right=.96,bottom=.15,top=.95,hspace=.25,wspace=.3)

    savefig("figs/4.eps")
end

# ----------------
# Figure 5
# ----------------
if in(5,plot_these)
    figure(5,(15,8))

    Ts = 282:.2:307
    λ = -.88
    N(T,F,b=-.0001,g=0,a=.058) = F + λ * (T - 287) + a * (T - 287)^2 + b * (T - 287)^3 + g * (T - 287)^5

    function model_setup(F,i)
        annotate(string("abcdef"[i]),xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
        g = -4e-6
        plot(Ts, [N(T,F,0,0,0) for T=Ts], linestyle = "--", color="0.5", linewidth=2)
        plot(Ts, [N(T,F,0,0) for T=Ts], color="0.5", linewidth=2)
        plot(Ts, [N(T,F,0,g) for T=Ts], linewidth=2)
        ylim(-4,4)
        xlim(Ts[1],Ts[end])
        # x axis
        plot(Ts,[0 for T=Ts],"k",linewidth=2)
        plot([287,287],[-.2,.2],"k",linewidth=2)
        xlabel("\$T\$ (\$K\$)",labelpad=-2.5)
        ylabel("\$N\$ (\$W/m^2\$)",labelpad=-10)
    
        ms = 8
        if F > 0
            arrow(287,.4,0,F - 1,shape="full",facecolor="k",head_width=.4,head_length=.2)
            text(285.3,F/2 - .3,"\$F\$")
        
            ΔT_is = ΔT(F,λ,.058,g)
            arrow(287.7,-.2,ΔT_is-1.6,0,shape="full",edgecolor="k",facecolor="k",head_width=.2,head_length=.4)
            text(287 + ΔT_is/2 - 1.2,-.9,"\$\\Delta T\$")
            if F == 3.28
               ΔT_2 = sort(fzeros(ΔT -> F + λ * ΔT + .058 * ΔT^2 + g * ΔT^5, 0, 50))[2] 
               plot(ΔT_2 + 287,0,marker="o",fillstyle="none",color="k",markersize=ms,mew=1.5)
           
               ΔT_3 = sort(fzeros(ΔT -> F + λ * ΔT + .058 * ΔT^2 + g * ΔT^5, 0, 50))[3]
               plot(ΔT_3 + 287,0,marker="o",fillstyle="full",color="k",markersize=ms)
            end
            if F == 3.46
                plot(ΔT_is + 287,0,marker="o",fillstyle="left",color="k",markersize=ms,mew=1.5)
                ΔT_2 = sort(fzeros(ΔT -> F + λ * ΔT + .058 * ΔT^2 + g * ΔT^5, 0, 50))[3] 
                plot(ΔT_2 + 287,0,marker="o",fillstyle="full",color="k",markersize=ms,mew=1.5)
            else
                plot(ΔT_is + 287,0,marker="o",fillstyle="full",color="k",markersize=ms)
                if F == 3.105
                    ΔT_2 = sort(fzeros(ΔT -> F + λ * ΔT + .058 * ΔT^2 + g * ΔT^5, 0, 50))[2] 
                    plot(ΔT_2 + 287,0,marker="o",fillstyle="right",color="k",markersize=ms,mew=1.5) 
                end
            end
        
            text(Ts[1]+.5,-3.8,"\$F = $(round(F,1)) W/m^2  \\Delta T = $(round(ΔT_is,1)) K\$")
        else
            text(286,-1,"\$T_0\$")
            plot(287,0,marker="o",fillstyle="full",color="k",markersize=ms)        
        end
    end

    subplot(231)
    model_setup(0,1)

    subplot(232)
    model_setup(2,2)

    subplot(233)
    model_setup(3.105,3)

    subplot(234)
    model_setup(3.28,4)

    subplot(235)
    model_setup(3.46,5)

    subplot(236)
    model_setup(3.6,6)

    subplots_adjust(left=.06,right=.96,bottom=.1,top=.95,hspace=.25)
    savefig("figs/5.eps")
end


# ----------------
# Figure 6
# ----------------
if in(6,plot_these)
    figure(6,(15,5))

    # GCM data
    models = [
        {
            :name => "Stouffer and Manabe\n(2003) - GFDL",
            :T_s => [-7.85, 0.0, 4.33, 7.9],
            :Net => [3.71,0.0,-3.71,-7.42],
            :D => [-1,0,1,2],
            :T_2x => 4.33,
            :T_4x => 7.9,
            :T_8x => 10000 # value unknown
        },
        {
            :name => "Hansen et al.\n(2005) - GISS",
            :T_s => [-4.56, -3.32, -1.78, .0, 0.58, 1.10, 1.96, 4.06, 7.02],
            :Net => [12.68,8.88,4.61,-0.0,-1.44,-2.64,-4.52,-9.27,-14.65],
            :D => [-3,-2,-1,0,.322,.585,1,2,3],
            :T_2x => 1.96,
            :T_4x => 4.06,
            :T_8x => 7.02
        },
        {
            :name => "Colman and McAvaney\n(2009) - BMRC",
            :T_s => [-12.3875, -7.74805, -4.60923, -2.70842, 0., 2.49182, 5.11104, 6.84392, 8.9713, 12.7421],
            :Net => [15.392,11.704,7.936,4.04,-0.0,-4.28,-8.8,-13.776,-19.848,-26.24],
            :D => [-4,-3,-2,-1,0,1,2,3,4,5],
            :T_2x => 2.49182,
            :T_4x => 5.11104,
            :T_8x => 6.84392
        },
        {
            :name => "Meraner et al.\n(2013) - ECHAM6",
            :T_s => [0,2.79,6.7,12.46,22.68],
            :Net => [0,-4.26,-9.09,-14.24,-19.18],
            :D => [0,1,2,3,4],
            :T_2x => 2.79,
            :T_4x => 6.7,
            :T_8x => 12.46
        }
    ]

    # helper functions/values
    function polyfit_exps_CO2(x,y,D)
        p = [1,1,2]
        A = float([ (j==2 ? D[i] : 1) * float(x[i])^p[j] for i = 1:length(x), j = 1:3 ])
        A \ y
    end

    function run_model(m,F)
        markersize = (m[:name][1:2] == "Me" ? 21 : (m[:name][1] == 'H' ? 17 : 15))
        λ, b = [round(x,3) for x=polyfit_exps_CO2(m[:T_s], m[:Net], m[:D])]
        # λ, a = [round(x,3) for x=polyfit_exps(m[:T_s], m[:Net],2,:no_constant)]
        plot(b,(F==3.71? m[:T_2x] : (F==7.42 ? m[:T_4x] : m[:T_8x]) ),marker="\$$(m[:name][1])\$",markersize=markersize,linestyle="none",color="k",label=m[:name])
    end

    # range of λs from CMIP5, assuming linearity (Andrews et al., 2012)
    CMIP5_range = {
        :max => [-.79, "purple"],
        :mean => [-1.17, "orange"],
        :min => [-1.78, "green"]
    }

    bs = -.15:.01:.15

    function subplot_setup(i)
        # helper values and functions
        doublings = 2^i
        F = 3.71i
        y_max = 10.5i
        function plot_pair(kv_pair,F)
            stat, pair = kv_pair
            λ, color = pair
            plot(0,-F/λ,marker=".",color=color,markersize=20,linestyle="none")
            plot(bs,[ΔT_w_CO2(i,λ,0,b) for b=bs], color=color, label="\$\\lambda_{$stat} \$ = $λ", linewidth=2)
        end

        # plot GCMs
        plot(1000,1000,linestyle="none",label="GCM experiments:")
        [run_model(m,F) for m=models]
        # plot lines
        plot(1000,1000,linestyle="none",label=" ")
        plot(1000,1000,linestyle="none",label="CMIP5 \$\\lambda\$s:")
        [plot_pair(pair,F) for pair=CMIP5_range]

        # center line
        plot([0 for f=[0,y_max]],[0,y_max],color="0.4")

        # labels and axes
        annotate(string("abc"[i]),xy=(subplot_label_x,subplot_label_y),xycoords="axes fraction",fontsize=16,horizontalalignment="left",verticalalignment="bottom")
        title("\$$(doublings)xCO_2\$")
        xlabel("\$b\$ (\$W/m^2/K\$ per doubling)")
        xlim(bs[1],bs[end])
        ylabel("\$\\Delta T_{$(doublings)x}\$ (\$K\$)",labelpad=-5)
        ylim(0,y_max)
    
        xticks(fontsize=14)
    end

    subplot(131)

        subplot_setup(1)

    subplot(132)

        subplot_setup(2)

    subplot(133)

        subplot_setup(3)
        legend(loc="center left", bbox_to_anchor=(1.1,.43),numpoints=1,prop={"size"=>13},frameon=false)


    subplots_adjust(left=.06,right=.8,bottom=.2,top=.85,wspace=.26)
    savefig("figs/6.eps")
end

