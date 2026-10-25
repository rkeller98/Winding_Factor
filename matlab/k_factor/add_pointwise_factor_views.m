function add_pointwise_factor_views(fig,x,y,groups,excitation,reduction,set_id)
% Pointwise alignment and plots; the conventional map/meshes are reused as-is.
map = reduction.map;
baseline = reduction.psi_reconstructed;
triangles = reduction.triangles;
tabs = findobj(fig,'Type','uitabgroup','Tag','contour_tabs');
tab = uitab(tabs,'Title','Pointwise k-factor','Tag','pointwise_factor');
note = uicontrol(tab,'Style','text','Units','normalized', ...
    'Position',[.03 .01 .77 .18],'HorizontalAlignment','left','FontSize',10, ...
    'FontName','Consolas','String','Open this tab to calculate fixed-map pointwise alignment.');
ready = false;
busy = false;
shown = true(numel(excitation),1);
checks = gobjects(numel(excitation),1);
views = gobjects(5,1);
built = false(5,1);
current = 1;
error_mode = false;
result = [];
mode_control = gobjects(0);
heat_control = gobjects(0);
dq_control = gobjects(0);
dq_mode = 1;
global_mode = false;
diagnostic_mode = false;
four_mode = false;
parameter_control = gobjects(0);
dq_busy = false;
heat_group = 1;
ie_limits = [min(x(:,3)),max(x(:,3))];
error_colors = [linspace(.1,1,128)' linspace(.3,1,128)' ones(128,1); ...
    ones(128,1) linspace(1,.2,128)' linspace(1,.1,128)'];
tabs.SelectionChangedFcn = @select_tab;

    function select_tab(~,event)
        if event.NewValue~=tab || ready || busy
            return
        end
        busy = true;
        note.String = 'Calculating bounded pointwise alignment on the fixed conventional map...';
        drawnow;
        try
            result = calculate_alignment(x,y,map);
        catch exception
            busy = false;
            note.String = ['Pointwise alignment unavailable: ' exception.message];
            warning('WindingFactor:PointwiseAlignment','%s',exception.message);
            return
        end
        % Numerical comparison uses identical physical points for both methods.
        result.comparison = result.valid_rec & all(isfinite([y baseline]),2);
        result.folded = cell(numel(excitation),1);
        result.folded_count = 0;
        for g = 1:numel(excitation)
            rows = find(groups==g);
            tri = triangles{g};
            supported = all(result.valid_rec(rows(tri)),2);
            before = signed_area(x(rows,1:2),tri);
            after = signed_area([result.shifted(rows),x(rows,2)],tri);
            % Orientation reversal or near-collapse of an ORIGINAL triangle.
            result.folded{g} = supported & after./before<=1e-8;
            result.folded_count = result.folded_count+nnz(result.folded{g});
        end
        ready = true;
        busy = false;
        controls = uipanel(tab,'Title','Excitation groups','Units','normalized', ...
            'Position',[.83 .02 .16 .96]);
        uicontrol(controls,'Style','pushbutton','String','All','Tag','pointwise_all', ...
            'Units','normalized','Position',[.04 .92 .43 .05], ...
            'Callback',@(~,~) show_all(true));
        uicontrol(controls,'Style','pushbutton','String','None','Tag','pointwise_none', ...
            'Units','normalized','Position',[.53 .92 .43 .05], ...
            'Callback',@(~,~) show_all(false));
        height = min(.05,.65/numel(excitation));
        for g = 1:numel(excitation)
            rows = groups==g;
            checks(g) = uicontrol(controls,'Style','checkbox','String', ...
                sprintf('%g A ref (%d/%d valid)',excitation(g),nnz(rows & result.valid_rec),nnz(rows)), ...
                'Value',1,'FontSize',9,'Units','normalized', ...
                'Position',[.04 .87-g*height .94 height], ...
                'Tag',sprintf('pointwise_check_%02d',g), ...
                'Callback',@(source,~) show_group(g,logical(source.Value)));
        end
        mode_control = uicontrol(controls,'Style','popupmenu','String', ...
            {'Measured + pointwise','Signed errors + heatmaps'},'Tag','pointwise_error_mode', ...
            'Units','normalized','Position',[.04 .10 .92 .05], ...
            'Visible','off','Callback',@change_mode);
        candidates = find(accumarray(groups,double(result.valid_rec))>0);
        if ~isempty(candidates)
            [~,nearest] = min(abs(excitation(candidates)-median(excitation)));
            heat_group = candidates(nearest);
        end
        heat_control = uicontrol(controls,'Style','popupmenu','String', ...
            compose('Heatmap: %g A ref',excitation),'Value',heat_group, ...
            'Tag','pointwise_heat_group','Units','normalized', ...
            'Position',[.04 .16 .92 .045],'Visible','off','Callback',@change_heat);
        dq_control = uicontrol(controls,'Style','popupmenu','String', ...
            {'Shifted measured maps','Fixed reduced reference map', ...
             'Measured + reconstruction','Signed errors + heatmaps'}, ...
            'Tag','twoaxis_view_mode','Units','normalized','Position',[.04 .10 .92 .05], ...
            'Visible','off','Callback',@change_dq_mode);
        parameter_control=uicontrol(controls,'Style','popupmenu', ...
            'String',{'Pointwise fixed-map','Consistent global','Fixed-map global diagnostic', ...
                      'Global 4-parameter transformation'}, ...
            'Tag','twoaxis_parameter_mode', ...
            'Units','normalized','Position',[.04 .05 .92 .045], ...
            'Visible','off','Callback',@change_parameter_mode);
        uicontrol(controls,'Style','text','String','Rotate / zoom with mouse.', ...
            'Units','normalized','Position',[.04 .005 .92 .035]);
        child = uitabgroup(tab,'Units','normalized','Position',[0 .21 .82 .79], ...
            'Tag','pointwise_views');
        names = {'Pointwise factor','Pointwise shifted maps','Pointwise reconstruction', ...
            'Pointwise k_d / k_q factors','Global / pointwise k_d / k_q reconstruction'};
        for k = 1:5
            views(k) = uitab(child,'Title',names{k},'UserData',k);
        end
        child.SelectionChangedFcn = @(~,event) draw_view(event.NewValue.UserData);
        draw_view(1);
        fprintf(['Set%d pointwise: %d/%d valid factors; %d common comparison points; ' ...
            '%d folded/collapsed shifted triangles marked.\n'],set_id,nnz(result.valid_rec), ...
            size(x,1),nnz(result.comparison),result.folded_count);
        fprintf('  RMSE constant/pointwise [d q joint] mVs: [%g %g %g] / [%g %g %g]\n', ...
            result.rmse_constant,result.rmse_pointwise);
        rotate3d(fig,'on');
    end

    function show_all(visible)
        shown(:) = visible;
        set(checks,'Value',double(visible));
        update_selection();
    end

    function show_group(g,visible)
        shown(g) = visible;
        update_selection();
    end

    function update_selection()
        for g = 1:numel(excitation)
            state = 'off';
            if shown(g)
                state = 'on';
            end
            set(findall(tab,'Tag',sprintf('pointwise_group_%02d',g)),'Visible',state);
        end
        update_note();
    end

    function change_mode(source,~)
        error_mode = source.Value==2;
        built(3) = false;
        draw_view(3);
    end

    function change_heat(source,~)
        heat_group = source.Value;
        built([3 5]) = false;
        if current>=4, draw_view(5); else, draw_view(3); end
    end

    function change_dq_mode(source,~)
        dq_mode = source.Value;
        built(5) = false;
        draw_view(5);
    end

    function change_parameter_mode(source,~)
        global_mode=source.Value>=2;
        diagnostic_mode=source.Value==3;
        four_mode=source.Value==4;
        built(5)=false;
        draw_view(5);
    end

    function draw_view(k)
        current = k;
        child=findobj(tab,'Type','uitabgroup','Tag','pointwise_views');
        note.Position=[.03 .01 .77 .18];
        child.Position=[0 .21 .82 .79];
        mode_control.Visible = 'off';
        heat_control.Visible = 'off';
        dq_control.Visible = 'off';
        parameter_control.Visible = 'off';
        if k>=4
            note.Position=[.03 .005 .77 .24];
            child.Position=[0 .26 .82 .74];
            if ~prepare_twoaxis(), return; end
            if k==5 && four_mode && ~prepare_four(), return; end
            if k==5
                dq_control.Visible = 'on';
                parameter_control.Visible = 'on';
                if dq_mode==4, heat_control.Visible = 'on'; end
            end
            if ~built(k)
                draw_twoaxis(k);
                built(k)=true;
            end
            update_selection();
            return
        end
        if k==3
            mode_control.Visible = 'on';
            if error_mode
                heat_control.Visible = 'on';
            end
        end
        if built(k)
            update_note();
            return
        end
        delete(views(k).Children);
        if k==1
            layout = tiledlayout(views(k),1,1,'TileSpacing','compact','Padding','compact');
            ax = nexttile(layout);
            hold(ax,'on');
            for g = 1:numel(excitation)
                rows = find(groups==g);
                tri = triangles{g};
                tri = tri(all(result.valid_rec(rows(tri)),2),:);
                z = result.k_rec(rows);
                tag = sprintf('pointwise_group_%02d',g);
                trisurf(tri,x(rows,1),x(rows,2),z,x(rows,3),'Parent',ax, ...
                    'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                rows = rows(result.valid_rec(rows));
                scatter3(ax,x(rows,1),x(rows,2),result.k_rec(rows),8,x(rows,3), ...
                    'filled','Tag',tag);
            end
            decorate(ax,'i_d in A','k^*_{rec} in A/A','Pointwise fixed-map alignment factor',false);
        else
            if k==3 && error_mode
                layout = tiledlayout(views(k),2,2,'TileSpacing','compact','Padding','compact');
            else
                layout = tiledlayout(views(k),1,2,'TileSpacing','compact','Padding','compact');
            end
            for c = 1:2
                if k==3 && error_mode
                    ax = nexttile(layout,(c-1)*2+1);
                    heat = nexttile(layout,(c-1)*2+2);
                    hold(heat,'on');
                else
                    ax = nexttile(layout,c);
                end
                hold(ax,'on');
                for g = 1:numel(excitation)
                    rows = find(groups==g);
                    original = triangles{g};
                    supported = all(result.valid_rec(rows(original)),2);
                    tri = original(supported,:);
                    tag = sprintf('pointwise_group_%02d',g);
                    coords = x(rows,1:2);
                    measured = 1000*y(rows,c);
                    predicted = 1000*result.prediction(rows,c);
                    if k==2
                        coords(:,1) = result.shifted(rows);
                        trisurf(tri,coords(:,1),coords(:,2),measured,x(rows,3),'Parent',ax, ...
                            'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                        folded = original(result.folded{g},:);
                        if ~isempty(folded)
                            trisurf(folded,coords(:,1),coords(:,2),measured,'Parent',ax, ...
                                'FaceColor','none','EdgeColor',[.8 0 .65], ...
                                'LineWidth',1.1,'LineStyle','--','Tag',tag);
                        end
                    elseif ~error_mode
                        % Keep the measured surface, including unsupported pointwise regions.
                        trisurf(original,coords(:,1),coords(:,2),measured,x(rows,3),'Parent',ax, ...
                            'FaceColor','interp','FaceAlpha',.5,'EdgeColor','none','Tag',tag);
                        trisurf(tri,coords(:,1),coords(:,2),predicted,'Parent',ax, ...
                            'FaceColor','none','EdgeColor',[.12 .12 .12], ...
                            'LineWidth',.45,'Tag',tag);
                    else
                        error = predicted-measured;
                        trisurf(tri,coords(:,1),coords(:,2),error,error,'Parent',ax, ...
                            'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                        if g==heat_group
                            trisurf(tri,coords(:,1),coords(:,2),zeros(size(error)),error,'Parent',heat, ...
                                'FaceColor','interp','EdgeColor','none','Tag',tag);
                            good = result.valid_rec(rows);
                            scatter(heat,coords(good,1),coords(good,2),8,error(good),'filled','Tag',tag);
                            scatter(heat,coords(~good,1),coords(~good,2),12,[.45 .45 .45], ...
                                'x','Tag',tag); % Unavailable locations, no fabricated error value.
                        end
                    end
                end
                symbol = sprintf('\\psi_%s',char('d'+13*(c-1)));
                if k==2
                    decorate(ax,'s_P = i_d + k^*_{rec} i_e in A',[symbol ' in mVs'], ...
                        [symbol ': shifted measured values'],false);
                elseif ~error_mode
                    decorate(ax,'i_d in A',[symbol ' in mVs'],[symbol ': measured + pointwise'],false);
                    h1 = plot3(ax,nan,nan,nan,'s','MarkerFaceColor',[.2 .5 .9],'MarkerEdgeColor','none');
                    h2 = plot3(ax,nan,nan,nan,'-','Color',[.12 .12 .12]);
                    legend(ax,[h1 h2],{'Measured (surface)','Pointwise (mesh)'}, ...
                        'Location','southoutside','FontSize',9,'AutoUpdate','off');
                else
                    decorate(ax,'i_d in A','Pointwise - measured in mVs',[symbol ': signed error'],true);
                    decorate(heat,'i_d in A','', ...
                        sprintf('%s: error at %g A ref',symbol,excitation(heat_group)),true);
                    view(heat,2);
                end
            end
        end
        built(k) = true;
        update_selection();
    end

    function ok = prepare_twoaxis()
        ok=isfield(result,'twoaxis');
        if ok || dq_busy, return; end
        dq_busy=true;
        note.String='Calculating global/local alignment on the fixed two-axis reference...';
        drawnow;
        try
            r=calculate_twoaxis_alignment(x,y,set_id,[map.p 0]);
        catch exception
            dq_busy=false;
            note.String=['Two-axis alignment unavailable: ' exception.message];
            warning('WindingFactor:TwoAxisAlignment','%s',exception.message);
            return
        end
        r.comparison=r.valid_rec & r.global.valid_rec & result.comparison;
        r.folded=cell(numel(excitation),1);
        r.folded_count=0;
        r.global.folded=cell(numel(excitation),1);
        r.global.folded_count=0;
        r.fixed_global.folded=cell(numel(excitation),1);
        r.fixed_global.folded_count=0;
        for g=1:numel(excitation)
            rows=find(groups==g);
            tri=triangles{g};
            supported=all(r.valid_rec(rows(tri)),2);
            before=signed_area(x(rows,1:2),tri);
            after=signed_area(r.shifted(rows,:),tri);
            r.folded{g}=supported & after./before<=1e-8;
            r.folded_count=r.folded_count+nnz(r.folded{g});
            supported=all(r.global.valid_rec(rows(tri)),2);
            after=signed_area(r.global.shifted(rows,:),tri);
            r.global.folded{g}=supported & after./before<=1e-8;
            r.global.folded_count=r.global.folded_count+nnz(r.global.folded{g});
            supported=all(r.fixed_global.valid_rec(rows(tri)),2);
            after=signed_area(r.fixed_global.shifted(rows,:),tri);
            r.fixed_global.folded{g}=supported & after./before<=1e-8;
            r.fixed_global.folded_count=r.fixed_global.folded_count+nnz(r.fixed_global.folded{g});
        end
        result.twoaxis=r;
        dq_busy=false;
        ok=true;
        common=r.comparison;
        summary=comparison_rmse(common);
        fprintf(['Set%d two-axis: reference %s [%.9g %.9g] A/A; %d/%d valid; ' ...
            'nonunique %d, weak excitation %d, ill-conditioned %d, unsupported %d.\n'], ...
            set_id,r.reference_name,r.map.p,nnz(r.valid_rec),size(x,1),nnz(r.nonunique), ...
            nnz(r.low_excitation),nnz(r.ill_conditioned),nnz(r.unsupported));
        fprintf('  Four-way common support: %d/%d. RMSE [d q joint] mVs:\n',nnz(common),size(x,1));
        fprintf('  Constant scalar: [%g %g %g]\n',summary(1,:));
        fprintf('  Consistent global kd/kq: [%g %g %g]\n',summary(2,:));
        fprintf('  Pointwise scalar: [%g %g %g]\n',summary(3,:));
        fprintf('  Pointwise kd/kq: [%g %g %g]\n',summary(4,:));
        e=1000*r.error(r.valid_rec,:);
        fprintf('  All valid two-axis points: RMSE [%g %g %g] mVs.\n', ...
            sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all')));
        g=r.global;
        pair_common=g.valid_rec & r.valid_rec;
        eg=1000*g.error(pair_common,:); ep=1000*r.error(pair_common,:);
        fprintf('  Consistent global fit [%.9g %.9g]; supported %d/%d (%.2f%%).\n', ...
            g.p,nnz(g.valid_rec),size(x,1),100*nnz(g.valid_rec)/size(x,1));
        fprintf('  Global/local pair common %d: RMSE global [%g %g %g], local [%g %g %g] mVs.\n', ...
            nnz(pair_common),sqrt(mean(eg.^2,1)),sqrt(mean(eg.^2,'all')), ...
            sqrt(mean(ep.^2,1)),sqrt(mean(ep.^2,'all')));
        fprintf(['  Consistent objective: fixed one-cell-interior N=%d; ANY lost fit point infeasible; no missing-data residual penalty.\n' ...
            '  Seed/fitted mean flux objective %.12g / %.12g Vs^2; %d lost population points.\n' ...
            '  Bounds kd=[%.6g %.6g], kq=[%.6g %.6g]; %d starts; parameter spread [%.6g %.6g], score spread %.9g.\n'], ...
            nnz(g.population),g.seed_objective,g.objective,nnz(g.population & ~g.valid_rec),g.bounds(:,1),g.bounds(:,2), ...
            size(g.start_pairs,1),g.initialization_spread,g.score_spread);
        for j=1:size(g.start_pairs,1)
            fprintf('    Start [%.8g %.8g] -> [%.8g %.8g], flux cost %.12g, flag %d, evaluations %d.\n', ...
                g.start_pairs(j,:),g.solutions(j,:),g.scores(j),g.exitflags(j),g.evaluations(j));
        end
        d=r.fixed_global;
        missing_term=2*nnz(d.population & ~d.valid_rec)/nnz(d.population);
        fprintf('  Preserved fixed-map diagnostic [%.9g %.9g]: cost %.12g, missing term %.12g (%.4f%%), initialization spread [%.6g %.6g].\n', ...
            d.p,d.objective,missing_term,100*missing_term/d.objective,d.initialization_spread);
        fprintf('  Pointwise geometric feasibility %d/%d; solver completed %d; failures %d; weak ie %d, multiple %d, ill-conditioned %d; uniquely retained %d.\n', ...
            nnz(r.geometric_supported),size(x,1),nnz(r.solver_completed),nnz(r.optimization_failed), ...
            nnz(r.low_excitation),nnz(r.nonunique),nnz(r.ill_conditioned),nnz(r.valid_rec));
        fprintf('  SAME-F_ref diagnostic supported but pointwise filtered: %d; infeasible among diagnostic-supported: %d.\n', ...
            nnz(d.valid_rec & ~r.valid_rec),nnz(d.valid_rec & ~r.geometric_supported));
        if ~isempty(r.historical_cv)
            h=r.historical_cv;
            fprintf('  ARCHIVED method-F excitation-group holdout: RMSE [%g %g %g] mVs, N=%d (%.2f%%). Reused old procedure; NOT validation of the refined global fit.\n', ...
                h.rmse,h.points,100*h.coverage);
        end
    end

    function summary = comparison_rmse(common)
        summary=nan(4,3);
        if ~any(common), return; end
        prediction={baseline,result.twoaxis.global.prediction,result.prediction,result.twoaxis.prediction};
        for method=1:4
            e=1000*(prediction{method}(common,:)-y(common,:));
            summary(method,:)=[sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all'))];
        end
    end

    function ok = prepare_four()
        ok=isfield(result.twoaxis,'four');
        if ok, return; end
        note.String='Fitting the global four-parameter transformation and its matching LUT...';
        drawnow;
        try
            result.twoaxis.four=fit_global_fourparameter(x,y,groups,triangles,result.twoaxis.global);
            built(5)=false;
            ok=true;
        catch exception
            note.String=['Four-parameter transformation unavailable: ' exception.message];
            warning('WindingFactor:FourParameter','%s',exception.message);
        end
    end

    function draw_twoaxis(k)
        delete(views(k).Children);
        r=result.twoaxis;
        label='Pointwise';
        if k==5 && global_mode
            r=result.twoaxis.global;
            label='Consistent global';
            if diagnostic_mode
                r=result.twoaxis.fixed_global;
                label='Fixed-map diagnostic';
            elseif four_mode
                r=result.twoaxis.four;
                label='Global 4-parameter';
            end
        end
        if k==5 && dq_mode==4
            layout=tiledlayout(views(k),2,2,'TileSpacing','compact','Padding','compact');
        else
            layout=tiledlayout(views(k),1,2,'TileSpacing','compact','Padding','compact');
        end
        error_limit=max(1000*abs([result.twoaxis.error;result.twoaxis.global.error;result.twoaxis.fixed_global.error]),[],'all','omitnan');
        if isfield(result.twoaxis,'four')
            error_limit=max(error_limit,max(1000*abs(result.twoaxis.four.error),[],'all','omitnan'));
        end
        for c=1:2
            if k==5 && dq_mode==4
                ax=nexttile(layout,(c-1)*2+1);
                heat=nexttile(layout,(c-1)*2+2);
                hold(heat,'on');
            else
                ax=nexttile(layout,c);
            end
            hold(ax,'on');
            symbol=sprintf('\\psi_%s',char('d'+13*(c-1)));
            if k==5 && dq_mode==2
                sd=r.map.origin(1)+r.map.scale(1)*r.map.axes{1};
                sq=r.map.origin(2)+r.map.scale(2)*r.map.axes{2};
                [D,Q]=ndgrid(sd,sq);
                surf(ax,D,Q,1000*r.map.values{c},'FaceColor','none','EdgeColor','interp','LineWidth',.5);
                xlabel(ax,'s_d in A'); ylabel(ax,'s_q in A'); zlabel(ax,[symbol ' in mVs']);
                if global_mode && ~diagnostic_mode
                    title(ax,[symbol ': matching ' label ' LUT'],'FontSize',12);
                else
                    title(ax,[symbol ': fixed pointwise reference'],'FontSize',12);
                end
                grid(ax,'on'); view(ax,35,25); colormap(ax,parula(256));
                both=1000*[result.twoaxis.map.values{c}(:);result.twoaxis.global.map.values{c}(:)];
                if isfield(result.twoaxis,'four')
                    both=[both;1000*result.twoaxis.four.map.values{c}(:)];
                end
                clim(ax,[min(both,[],'omitnan'),max(both,[],'omitnan')]);
                cb=colorbar(ax); cb.Label.String=[symbol ' in mVs'];
                if global_mode && ~diagnostic_mode && isfield(result.twoaxis,'four')
                    coords=[result.twoaxis.global.shifted;result.twoaxis.four.shifted];
                    xlim(ax,[min(coords(:,1)) max(coords(:,1))]);
                    ylim(ax,[min(coords(:,2)) max(coords(:,2))]);
                    zlim(ax,[min(both,[],'omitnan') max(both,[],'omitnan')]);
                end
                continue
            end
            for g=1:numel(excitation)
                rows=find(groups==g);
                original=triangles{g};
                tri=original(all(r.valid_rec(rows(original)),2),:);
                tag=sprintf('pointwise_group_%02d',g);
                coords=x(rows,1:2);
                measured=1000*y(rows,c);
                predicted=1000*r.prediction(rows,c);
                if k==4
                    z=r.factors(rows,c);
                    trisurf(tri,coords(:,1),coords(:,2),z,x(rows,3),'Parent',ax, ...
                        'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                    good=rows(r.valid_rec(rows));
                    scatter3(ax,x(good,1),x(good,2),r.factors(good,c),8,x(good,3),'filled','Tag',tag);
                elseif dq_mode==1
                    coords=r.shifted(rows,:);
                    trisurf(tri,coords(:,1),coords(:,2),measured,x(rows,3),'Parent',ax, ...
                        'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                    folded=original(r.folded{g},:);
                    if ~isempty(folded)
                        trisurf(folded,coords(:,1),coords(:,2),measured,'Parent',ax, ...
                            'FaceColor','none','EdgeColor',[.8 0 .65],'LineWidth',1.1,'LineStyle','--','Tag',tag);
                    end
                elseif dq_mode==3
                    trisurf(original,coords(:,1),coords(:,2),measured,x(rows,3),'Parent',ax, ...
                        'FaceColor','interp','FaceAlpha',.5,'EdgeColor','none','Tag',tag);
                    trisurf(tri,coords(:,1),coords(:,2),predicted,'Parent',ax, ...
                        'FaceColor','none','EdgeColor',[.12 .12 .12],'LineWidth',.45,'Tag',tag);
                else
                    error=predicted-measured;
                    trisurf(tri,coords(:,1),coords(:,2),error,error,'Parent',ax, ...
                        'FaceColor','none','EdgeColor','interp','LineWidth',.5,'Tag',tag);
                    if g==heat_group
                        trisurf(tri,coords(:,1),coords(:,2),zeros(size(error)),error,'Parent',heat, ...
                            'FaceColor','interp','EdgeColor','none','Tag',tag);
                        good=r.valid_rec(rows);
                        scatter(heat,coords(good,1),coords(good,2),8,error(good),'filled','Tag',tag);
                        scatter(heat,coords(~good,1),coords(~good,2),12,[.45 .45 .45],'x','Tag',tag);
                    end
                end
            end
            if k==4
                name=char('d'+13*(c-1));
                decorate(ax,'i_d in A',sprintf('k_%s^* in A/A',name), ...
                    sprintf('Pointwise k_%s: fixed-map alignment',name),false);
                [D,Q]=meshgrid([min(x(:,1)),max(x(:,1))],[min(x(:,2)),max(x(:,2))]);
                plane=surf(ax,D,Q,repmat(result.twoaxis.global.p(c),2,2), ...
                    'FaceColor',[.3 .3 .3],'FaceAlpha',.18,'EdgeColor',[.15 .15 .15], ...
                    'LineStyle','--','Tag','global_factor_plane');
                legend(ax,plane,sprintf('Consistent global k_%s = %.6g A/A',name,result.twoaxis.global.p(c)), ...
                    'Location','southoutside','FontSize',9,'AutoUpdate','off');
            elseif dq_mode==1
                decorate(ax,'s_d = i_d + k_d i_e in A',[symbol ' in mVs'], ...
                    [symbol ': ' label ' shifted measured values'],false);
                ylabel(ax,'s_q = i_q + k_q i_e in A');
                if four_mode
                    xlabel(ax,'s_d = i_d + k_{dq} i_q + k_{de} i_e in A');
                    ylabel(ax,'s_q = k_{qd} i_d + i_q + k_{qe} i_e in A');
                end
            elseif dq_mode==3
                decorate(ax,'i_d in A',[symbol ' in mVs'],[symbol ': measured + ' label],false);
                h1=plot3(ax,nan,nan,nan,'s','MarkerFaceColor',[.2 .5 .9],'MarkerEdgeColor','none');
                h2=plot3(ax,nan,nan,nan,'-','Color',[.12 .12 .12]);
                legend(ax,[h1 h2],{'Measured (surface)',[label ' (mesh)']}, ...
                    'Location','southoutside','FontSize',9,'AutoUpdate','off');
            else
                decorate(ax,'i_d in A','Reconstructed - measured in mVs', ...
                    [symbol ': ' label ' signed error'],true,error_limit);
                decorate(heat,'i_d in A','',sprintf('%s: error at %g A ref',symbol,excitation(heat_group)),true,error_limit);
                view(heat,2);
                for error_ax=[ax heat]
                    cb=colorbar(error_ax);
                    cb.Label.String='Reconstructed - measured in mVs';
                end
            end
            if k==5 && global_mode && ~diagnostic_mode && isfield(result.twoaxis,'four')
                if dq_mode==1
                    coords=[result.twoaxis.global.shifted;result.twoaxis.four.shifted];
                else
                    coords=x(:,1:2);
                end
                xlim(ax,[min(coords(:,1)) max(coords(:,1))]);
                ylim(ax,[min(coords(:,2)) max(coords(:,2))]);
                if dq_mode==4
                    xlim(heat,[min(x(:,1)) max(x(:,1))]);
                    ylim(heat,[min(x(:,2)) max(x(:,2))]);
                    zlim(ax,[-error_limit error_limit]);
                else
                    flux_limits=1000*[y(:,c);result.twoaxis.global.prediction(:,c);result.twoaxis.four.prediction(:,c)];
                    zlim(ax,[min(flux_limits,[],'omitnan') max(flux_limits,[],'omitnan')]);
                end
            end
        end
    end

    function decorate(ax,x_label,z_label,name,is_error,error_limit)
        xlabel(ax,x_label);
        ylabel(ax,'i_q in A');
        zlabel(ax,z_label);
        title(ax,name,'FontSize',12);
        grid(ax,'on');
        view(ax,35,25);
        cb = colorbar(ax);
        if is_error
            if nargin<6
                limit = max(1000*abs(result.prediction-y),[],'all','omitnan');
            else
                limit = error_limit;
            end
            if ~isfinite(limit) || limit==0
                limit = 1;
            end
            clim(ax,[-limit limit]);
            colormap(ax,error_colors);
            cb.Label.String = 'Pointwise - measured in mVs';
        else
            clim(ax,ie_limits);
            colormap(ax,turbo(256));
            cb.Label.String = 'Measured excitation current i_e in A';
        end
    end

    function update_note()
        active_valid=result.valid_rec;
        if current>=4 && isfield(result,'twoaxis')
            active_valid=result.twoaxis.valid_rec;
            if current==5 && global_mode
                active_valid=result.twoaxis.global.valid_rec;
                if diagnostic_mode, active_valid=result.twoaxis.fixed_global.valid_rec; end
                if four_mode && isfield(result.twoaxis,'four'), active_valid=result.twoaxis.four.valid_rec; end
            end
        end
        for g=1:numel(excitation)
            rows=groups==g;
            checks(g).String=sprintf('%g A ref (%d/%d valid)',excitation(g), ...
                nnz(rows & active_valid),nnz(rows));
        end
        if current>=4 && isfield(result,'twoaxis')
            selected=shown(groups);
            if current==5 && four_mode && isfield(result.twoaxis,'four')
                f=result.twoaxis.four; g=result.twoaxis.global;
                common=selected & f.common;
                eg=1000*g.error(common,:); ef=1000*f.error(common,:);
                geometry=f.geometry;
                note.String=sprintf(['Set%d | Global 4-parameter: [kdq kde kqd kqe]=[%.6g %.6g %.6g %.6g]\n' ...
                    'Effective [kd kq]=[%.6g %.6g] A/A; det(A)=%.5g; cond(A)=%.5g; matching 81 x 81 LUT\n' ...
                    'Selected support two/four: %d/%d, %d/%d; identical comparison points: %d\n' ...
                    'In-sample RMSE in mVs          d             q             joint\n' ...
                    'Global two-factor       %11.4f   %11.4f   %11.4f\n' ...
                    'Global four-parameter   %11.4f   %11.4f   %11.4f\n' ...
                    'All groups: reversed=%d, collapsed=%d, self-overlap pairs=%d; original group meshes.\n' ...
                    'Fixed fit N=%d; no lost fit points. Existing folds retained; no new self-overlaps accepted.\n' ...
                    'Invertible shear adds no continuous-model freedom. Grid/interpolation effects; no physical uniqueness claim.'], ...
                    set_id,f.p,f.effective,geometry.determinant,geometry.condition, ...
                    nnz(selected & g.valid_rec),nnz(selected),nnz(selected & f.valid_rec),nnz(selected),nnz(common), ...
                    sqrt(mean(eg.^2,1)),sqrt(mean(eg.^2,'all')),sqrt(mean(ef.^2,1)),sqrt(mean(ef.^2,'all')), ...
                    geometry.reversed,geometry.collapsed,geometry.adjacent+geometry.nonadjacent,nnz(f.population));
                tab.UserData=result;
                return
            end
            common=selected & result.twoaxis.comparison;
            summary=comparison_rmse(common);
            result.twoaxis.rmse_common=summary;
            result.twoaxis.comparison_points=nnz(common);
            tab.UserData=result;
            detail='Factors with multiple inverse solutions or poorly conditioned Jacobians remain undefined.';
            if current==5 && dq_mode==1
                folded_count=result.twoaxis.folded_count;
                if global_mode, folded_count=result.twoaxis.global.folded_count; end
                if diagnostic_mode && global_mode, folded_count=result.twoaxis.fixed_global.folded_count; end
                detail=sprintf('Magenta dashed: %d reversed/near-collapsed shifted triangles; original connectivity.', ...
                    folded_count);
            elseif current==5 && dq_mode==2
                detail='Global uses its matching LUT; pointwise/diagnostic use the unchanged fixed reference. No rebuild on switching.';
            elseif current==5 && dq_mode==4
                detail=sprintf('Heatmaps: %g A ref; gray crosses: unavailable errors. Each method uses its labelled LUT.',excitation(heat_group));
            end
            if current==5 && global_mode && diagnostic_mode
                d=result.twoaxis.fixed_global;
                detail=sprintf('DIAGNOSTIC VIEW ONLY: fixed-map pair [%.5g %.5g]; original 1-Vs penalty retained. Table uses consistent global.',d.p);
            end
            g=result.twoaxis.global;
            note.String=sprintf(['Set%d | Consistent global [%.5g %.5g]; fixed pointwise reference %s [%.5g %.5g] A/A\n' ...
                'Global support %d/%d; pointwise feasible %d, identifiable %d/%d; four-way common %d/%d (%.1f%%)\n' ...
                'RMSE in mVs                 d             q             joint\n' ...
                'Constant scalar       %11.4f   %11.4f   %11.4f\n' ...
                'Consistent global     %11.4f   %11.4f   %11.4f\n' ...
                'Pointwise scalar      %11.4f   %11.4f   %11.4f\n' ...
                'Pointwise k_d/k_q     %11.4f   %11.4f   %11.4f\n' ...
                '%s\nGlobal: fixed fit population, matching F_g. Scalar: F_A. Local: F_ref and target flux. All main errors are in-sample.'], ...
                set_id,g.p,result.twoaxis.reference_name,result.twoaxis.map.p, ...
                nnz(selected & g.valid_rec),nnz(selected),nnz(selected & result.twoaxis.geometric_supported), ...
                nnz(selected & result.twoaxis.valid_rec),nnz(selected), ...
                nnz(common),nnz(selected),100*nnz(common)/max(nnz(selected),1), ...
                summary(1,:),summary(2,:),summary(3,:),summary(4,:),detail);
            return
        end
        selected = shown(groups);
        common = selected & result.comparison;
        total = nnz(selected);
        retained = nnz(common);
        result.rmse_constant = nan(1,3);
        result.rmse_pointwise = nan(1,3);
        if retained>0
            a = 1000*(baseline(common,:)-y(common,:));
            p = 1000*(result.prediction(common,:)-y(common,:));
            result.rmse_constant = [sqrt(mean(a.^2,1)),sqrt(mean(a.^2,'all'))];
            result.rmse_pointwise = [sqrt(mean(p.^2,1)),sqrt(mean(p.^2,'all'))];
        end
        result.comparison_points = retained;
        tab.UserData = result; % GUI state only; no pointwise result files.
        detail = 'Undefined factors remain gaps; no replacement values or new fitted map.';
        if current==2
            detail = sprintf(['Magenta dashed: %d reversed/near-collapsed original triangles (all groups). ' ...
                'Invalid factors have no s_P and are omitted.'],result.folded_count);
        elseif current==3 && error_mode
            detail = sprintf(['Heatmaps: %g A ref; gray crosses: unavailable errors. ' ...
                'Original triangle connectivity; the group checkbox controls visibility.'],excitation(heat_group));
        end
        note.String = sprintf(['Set%d | Selected: %d/%d valid factors | Same-point RMSE support: %d/%d (%.1f%%)\n' ...
            'RMSE in mVs              d             q             joint\n' ...
            'Constant k_A       %11.4f   %11.4f   %11.4f\n' ...
            'Pointwise          %11.4f   %11.4f   %11.4f   (joint = sqrt(mean((ed^2+eq^2)/2)))\n' ...
            '%s\nIn-sample alignment diagnostic; not independent prediction accuracy or an identified physical parameter.\n' ...
            'k* uses the measured target flux and is not automatically a factor map for unseen operating points.'], ...
            set_id,nnz(selected & result.valid_rec),total,retained,total,100*retained/max(total,1), ...
            result.rmse_constant,result.rmse_pointwise,detail);
    end
end

function area = signed_area(x,tri)
a = x(tri(:,2),:)-x(tri(:,1),:);
b = x(tri(:,3),:)-x(tri(:,1),:);
area = a(:,1).*b(:,2)-a(:,2).*b(:,1);
end

function r = calculate_alignment(x,y,map)
% At fixed iq each supported grid cell is a straight segment in (Psi_d,Psi_q).
% Project the measurement onto EVERY segment and select the global minimum.
n = size(x,1);
s = map.origin(1)+map.scale(1)*map.axes{1};
intervals = numel(s)-1;
ds = diff(map.axes{1});
% Sample inside each cell, avoiding the ambiguous NaN-neighbour behaviour
% exactly at knots. Two interior samples determine the entire linear segment.
interior = [map.axes{1}(1:end-1)+ds/3,map.axes{1}(1:end-1)+2*ds/3];
S = repmat(interior,n,1);
Q = repmat((x(:,2)-map.origin(2))/map.scale(2),1,2*intervals);
profile = nan(n,2*intervals,2);
F = cell(1,2);
for c = 1:2
    F{c} = griddedInterpolant(map.axes,map.values{c},'linear','none');
    profile(:,:,c) = reshape(F{c}(S(:),Q(:)),n,[]);
end
first = profile(:,1:intervals,:);
second = profile(:,intervals+1:end,:);
left = 2*first-second;
delta = 3*(second-first);
target = reshape(y,n,1,2);
denominator = sum(delta.^2,3);
u = sum((target-left).*delta,3)./denominator;
u(denominator==0) = 0; % A flat segment is screened for ambiguity below.
u = min(max(u,0),1);
prediction = left+u.*delta;
loss = sum((prediction-target).^2,3);
loss(~all(isfinite(left) & isfinite(delta),3)) = inf;
[best,index] = min(loss,[],2);
row = (1:n)';
selected = sub2ind(size(loss),row,index);
candidate_s = s(1:end-1)+u.*diff(s);
best_s = candidate_s(selected);
r.unsupported = ~isfinite(best);
% Equal minima at the SAME knot are one solution. Separated minimizers or a
% minimizing flat segment are ambiguous; never select one arbitrarily.
loss_tolerance = 128*eps(max(sum(y.^2,2),1e-12));
s_tolerance = 1e-8*(max(s)-min(s));
tied = isfinite(loss) & abs(loss-best)<=loss_tolerance;
flat = denominator <= loss_tolerance;
r.nonunique = any(tied & (abs(candidate_s-best_s)>s_tolerance | flat),2);
grid_slope = sqrt(denominator)./diff(s);
sensitivity = abs(x(:,3)).*grid_slope(selected);
r.excitation_threshold = .1*median(abs(x(:,3)));
r.low_excitation = abs(x(:,3))<=r.excitation_threshold;
usable = ~r.unsupported & ~r.low_excitation & isfinite(sensitivity);
r.grid_threshold = .1*median(sensitivity(usable));
r.low_sensitivity = ~isfinite(sensitivity) | sensitivity<=r.grid_threshold;
r.valid_rec = ~r.unsupported & ~r.nonunique & ~r.low_excitation & ~r.low_sensitivity;
r.k_rec = (best_s-x(:,1))./x(:,3);
% Evaluate the SAME s_P used by the shifted view directly. Undoing k_A via a
% synthetic d-current adds cancellation and can move a valid knot into NaNs.
r.shifted = x(:,1)+r.k_rec.*x(:,3);
q = (x(:,2)-map.origin(2))/map.scale(2);
aligned = nan(n,2);
for c = 1:2
    aligned(:,c) = F{c}((r.shifted-map.origin(1))/map.scale(1),q);
end
% At a support boundary, roundoff can still place s_P just outside the cell.
% Move only those failed evaluations a few ulps toward the winning interval,
% accepting the correction only within that interval and objective roundoff.
repair = r.valid_rec & ~all(isfinite(aligned),2);
if any(repair)
    lower = reshape(s(index),[],1);
    upper = reshape(s(index+1),[],1);
    step = 8*eps(max([abs(x(:,1)),abs(best_s),ones(n,1)],[],2))./abs(x(:,3)) ...
        +8*eps(max(abs(r.k_rec),1));
    trial_k = r.k_rec+sign(((lower+upper)/2-best_s)./x(:,3)).*step;
    trial_s = x(:,1)+trial_k.*x(:,3);
    trial = nan(n,2);
    for c = 1:2
        trial(:,c) = F{c}((trial_s-map.origin(1))/map.scale(1),q);
    end
    accept = repair & trial_s>=lower & trial_s<=upper & all(isfinite(trial),2) ...
        & abs(sum((trial-y).^2,2)-best)<=loss_tolerance;
    r.k_rec(accept) = trial_k(accept);
    r.shifted(accept) = trial_s(accept);
    aligned(accept,:) = trial(accept,:);
end
r.unsupported = r.unsupported | ~all(isfinite(aligned),2);
r.valid_rec = r.valid_rec & ~r.unsupported & isfinite(r.k_rec);
r.k_rec(~r.valid_rec) = nan;
r.error_mVs = 1000*vecnorm(aligned-y,2,2);
r.error_mVs(r.unsupported) = nan;
r.prediction = aligned;
r.prediction(~r.valid_rec,:) = nan;
r.shifted(~r.valid_rec) = nan;
end
