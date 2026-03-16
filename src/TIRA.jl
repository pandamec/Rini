module TIRA
    
    using CSV
    using DataFrames
    using Polynomials
    using CairoMakie
    using Statistics

    function import_TIRA(file_path)
        
        df=CSV.read(file_path,DataFrame) 
        df.Zeit = parse.(Float64, replace.(df.Zeit, "," => "."))
        df.Länge = parse.(Float64, replace.(df.Länge, "," => "."))
        df.Weg = parse.(Float64, replace.(df.Weg, "," => "."))
        df.Kraft = parse.(Float64, replace.(df.Kraft, "," => "."))
        df.dL_ORG = parse.(Float64, replace.(df.dL_ORG, "," => "."))

        return df

    end

    function filter_range(df::DataFrame, column::String, low, high)
        return filter(row -> low ≤ row[column] ≤ high, df)
    end

    function get_corresponding(df::DataFrame, search_col::Symbol, search_val, return_col::Symbol)
        row = findfirst(df[!, search_col] .== search_val)
        return row === nothing ? missing : df[row, return_col]
    end


    function setZero(BaseName,name,i,ϵ_zero)

        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))

        df_zero=filter_range(df,"dL_ORG",ϵ_zero[1],ϵ_zero[2])

        df_zero[!,:dL_ORG]=  df_zero[!,:dL_ORG] .- df_zero[1,:dL_ORG]
        df_zero[!,:Kraft]=df_zero[!,:Kraft].-df_zero[1,:Kraft]
        
        return df_zero
    end


    function computeE(X,Y,T)


        p_fit=fit(X,Y,1)
        y_fit = p_fit.(X)

        E=p_fit[1] #MPa
        residuals = y_fit .- Y   # residuals
        n = length(X)
        σ = sqrt(sum(residuals.^2) / (n - 2))
        
        #σ_res = std(residuals, corrected=true) 
        X̄ = mean(X)
        u_E = σ / sqrt(sum((X .- X̄).^2))
        e=(residuals./Y)*100
        e_max=maximum(e)
        u_E=u_E*100/E

        fit_X_T    = fit(T,X,1)
        strainRate = fit_X_T[1]*60*100

        df_sum= DataFrame(Zeit=T,Strain=X*100, StressFit=y_fit*1000) # s Percentage MPa 
        
        
        return E, u_E, e_max, df_sum, strainRate

    end

    function computeStressStrainRaw(BaseName,name,i,Temp,As)

        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=(df[!,:dL_ORG]) #Computed with the videoextensometer value
        df[!,:Temperature]=fill(Temp,nrow(df))
        
        return df
    end


    
    function computeStressStrain(df,As,Temp)


        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=(df[!,:dL_ORG]) #Computed with the videoextensometer value
        df[!,:Temperature]=fill(Temp,nrow(df))

        return df
    end


    function computeProperties(df,LinearRange,StressRange)

        T=df[!,:Zeit] #s
        X=df[!,:dL_ORG] #%
        Y=df[!,:Stress] #MPa
        df_linear=filter_range(df,"Stress",StressRange[1],StressRange[2])
        df_linear=filter_range(df_linear,"dL_ORG",LinearRange[1],LinearRange[2])
           
        T_linear=df_linear[!,:Zeit] #s
        X_linear=df_linear[!,:dL_ORG] # Percentage
        Y_linear=df_linear[!,:Stress] #MPa

        E_modul,u_E,e_max,df_fit,strainRate  = computeE(X_linear/100,Y_linear/1000,T_linear) #Strain, #MGPa
        sigmaB      =   maximum(Y)
        sigmaF      =   maximum(Y_linear)
        δ_B         =   get_corresponding(df,:Stress,maximum(df[!,:Stress]),:dL_ORG)

        sum=Dict("Probe"=>df[1,:Name],
                 " Strain Rate[%/min]"=>strainRate,
                 "E[GPa] "=>E_modul,
                 "u_E[%] "=>u_E,
                 "u_E[GPa] "=>u_E*E_modul/100,
                 "e_max[%]" =>e_max,
                 "sigma0_25% "=> sigmaF,
                 "sigmaMax"=> sigmaB,
                 "δ_B"=> δ_B,
                 "δ_F/δ_B"=>df_linear[end,:dL_ORG]*100/δ_B,
                 "Temperature [C]"=> df_linear[1,:Temperature],
                 "LinearRange [%]"=>[df_linear[1,:dL_ORG],df_linear[end,:dL_ORG]])

       return sum,df_fit

    end

    export import_TIRA,computeStressStrain,computeProperties, filter_range

end