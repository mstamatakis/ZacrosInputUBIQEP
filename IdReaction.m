function type = IdReaction(GasReact,SurfReact,GasProd,SurfProd)

type = 0;

if ~isempty(GasReact)
    if ~isempty(SurfReact)
        disp('Reaction unknown');
    else
        if ~isempty(GasProd);
            disp('Reaction unknown');
        else
            if ~isempty(SurfProd)
                if length(SurfProd)==1
                    type=1;
                end
                if length(SurfProd)==2
                    if strcmp(SurfProd(1),SurfProd(2))==1
                        type=2;
                    else
                        type=4;
                    end
                end
            else
                disp('No Prodocuts insert');
            end
        end
    end
else
    if ~isempty(SurfReact)
        if ~isempty(GasProd)
            disp('Reaction unknown');
        else
            if ~isempty(SurfProd)
                if length(SurfReact)==1
                    type=5;
                end
                if length(SurfReact)==2
                    type=6;
                end
            else
                disp('No Prodocuts insert');
            end
        end
    else
        disp ('No Reactants or Products insert for this reaction');
    end
end

end