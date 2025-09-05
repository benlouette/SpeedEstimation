function db = buildIMSDatabase(source, outDir)
% buildIMSDatabase - Télécharge/ingère l'IMS Bearing Data et crée une base simple.
%
% USAGE:
%   % 1) Depuis une URL (ZIP):
%   db = buildIMSDatabase('https://exemple.com/IMS-Rexnord%20Bearing%20Data.zip', 'C:\IMS_DB');
%
%   % 2) Depuis un dossier local déjà extrait:
%   db = buildIMSDatabase('C:\path\to\IMS-Rexnord Bearing Data', 'C:\IMS_DB');
%
% RESULTAT:
%   - db.index : table avec (Set, FileName, Timestamp, Fs, N, NChannels, RelMatPath, AbsTxtPath, Duration_s)
%   - Sauvegarde des .mat individuels sous outDir\data\setX\*.mat
%   - Un fichier d’index sauvegardé en .mat et .csv

    arguments
        source (1,1) string
        outDir (1,1) string
    end

    Fs = 20000;          % d'après la doc IMS (20 kHz)
    N_expected = 20480;  % points par fichier (1 seconde)
    tsFormat = 'yyyy.MM.dd.HH.mm.ss';

    if ~exist(outDir, 'dir'); mkdir(outDir); end
    dataRoot = fullfile(outDir, 'data');
    if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

    % 1) Récupération des fichiers (ZIP ou dossier)
    if startsWith(lower(source), ["http://","https://"])
        zipFile = fullfile(outDir, 'ims_data.zip');
        fprintf('Téléchargement de l''archive depuis %s ...\n', source);
        websave(zipFile, source);
        fprintf('Extraction de l''archive...\n');
        extractedRoot = fullfile(outDir, 'extracted');
        unzip(zipFile, extractedRoot);
        rootPath = extractedRoot;
    else
        % source = dossier local (racine de l’archive extraite ou un sous-dossier “Set x”)
        rootPath = source;
    end

    % 2) Lister tous les fichiers ASCII (les jeux IMS sont souvent en sous-dossiers “1st_test”, “2nd_test”, etc.)
    fileList = listAllAsciiFiles(rootPath);
    if isempty(fileList)
        error('Aucun fichier texte trouvé sous: %s', rootPath);
    end
    fprintf('Trouvé %d fichiers.\n', numel(fileList));

    % 3) Boucle de parsing + sauvegarde
    rows = [];  % pour table
    rows(1).Set = [];  %#ok<AGROW> (préallocation dynamique)
    nOk = 0; nErr = 0;
    for i = 1:numel(fileList)
        txtPath = fileList(i);

        % Déterminer le “Set” depuis le chemin si possible (sinon 0)
        setNo = detectSetNumber(txtPath);

        % Extraire le timestamp à partir du nom de fichier
        [ts, tsOK] = parseTimestampFromName(txtPath, tsFormat);
        if ~tsOK
            warning('Timestamp non parsé: %s', txtPath);
            ts = NaT; %#ok<NASGU>
        end

        % Lire la matrice (robuste aux tabulations/espaces/quotes)
        try
            A = robustReadAsciiMatrix(txtPath);
            % A : (N, nCols), double
        catch ME
            warning('Lecture échouée (%s): %s', ME.identifier, txtPath);
            nErr = nErr + 1;
            continue
        end

        % Vérifications minimales
        [N, nCols] = size(A);
        if N ~= N_expected
            % On n’arrête pas—on garde l’info
            warning('Taille inattendue [%d x %d] (attendu %d x 4|8) fichier: %s', N, nCols, N_expected, txtPath);
        end
        if ~ismember(nCols, [4 8])
            warning('Nombre de colonnes inhabituel: %d (fichier: %s)', nCols, txtPath);
        end

        % Chemin de sortie .mat (par set)
        setDir = fullfile(dataRoot, sprintf('set%d', setNo));
        if ~exist(setDir, 'dir'); mkdir(setDir); end

        % Nom .mat = même nom de base que le fichier texte
        [~, baseName, ~] = fileparts(txtPath);
        matName = baseName + ".mat";
        matAbsPath = fullfile(setDir, matName);
        relMatPath = string(strrep(matAbsPath, outDir + filesep, ''));

        data  = A;            %#ok<NASGU>  % pour save
        FsLoc = Fs;           %#ok<NASGU>
        Ts    = ts;           %#ok<NASGU>
        NLoc  = N;            %#ok<NASGU>
        NCh   = nCols;        %#ok<NASGU>

        % Sauvegarde .mat (v7.3 si tu veux HDF5 pour très grands volumes)
        save(matAbsPath, 'data', 'FsLoc', 'Ts', 'NLoc', 'NCh', '-v7');

        % Ajouter une ligne à l’index
        nOk = nOk + 1;
        rows(nOk).Set         = setNo;           %#ok<AGROW>
        rows(nOk).FileName    = string(baseName);
        rows(nOk).Timestamp   = ts;
        rows(nOk).Fs          = Fs;
        rows(nOk).N           = N;
        rows(nOk).NChannels   = nCols;
        rows(nOk).RelMatPath  = relMatPath;      % chemin relatif au dossier outDir
        rows(nOk).AbsTxtPath  = string(txtPath);
        rows(nOk).Duration_s  = N / Fs;
    end

    % 4) Construire la table index
    db.index = struct2table(rows);
    db.root  = outDir;

    % 5) Sauvegarder l’index
    save(fullfile(outDir, 'ims_index.mat'), 'db');
    writetable(db.index, fullfile(outDir, 'ims_index.csv'));

    fprintf('Terminé. %d fichiers OK, %d erreurs.\n', nOk, nErr);
end

%% --- Helpers ---

function files = listAllAsciiFiles(rootPath)
    % Recherche tous les fichiers qui ressemblent aux snapshots IMS (pas .mat, .zip etc.)
    d = dir(fullfile(rootPath, '**', '*'));
    files = strings(0,1);
    for k = 1:numel(d)
        if ~d(k).isdir
                files(end+1,1) = string(fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
        end
    end
end

function setNo = detectSetNumber(p)
    p = string(p);
    lowerp = lower(p);
    % Heuristiques courantes : “1st_test”, “2nd_test”, “3rd_test”
    if contains(lowerp, "1st") || contains(lowerp, "set1")
        setNo = 1;
    elseif contains(lowerp, "2nd") || contains(lowerp, "set2")
        setNo = 2;
    elseif contains(lowerp, "3rd") || contains(lowerp, "set3")
        setNo = 3;
    else
        setNo = 0; % inconnu
    end
end

function [ts, ok] = parseTimestampFromName(p, fmt)
    % Cherche un motif yyyy.MM.dd.HH.mm.ss dans le nom de fichier
    [~, base, ext] = fileparts(p);
    base = base + ext;
    m = regexp(base, '\d{4}\.\d{2}\.\d{2}\.\d{2}\.\d{2}\.\d{2}', 'match', 'once');
    if ~isempty(m)
        try
            ts = datetime(m, 'InputFormat', fmt, 'TimeZone', 'UTC');
            ok = true;
            return
        catch
        end
    end
    ts = NaT; ok = false;
end

function A = robustReadAsciiMatrix(txtPath)
    % Lecture robuste de matrices ASCII (tabulations/espaces), nombres éventuellement entre guillemets.
    % Tente readmatrix en 'text', puis fallback textscan/sscanf si besoin.

    try
        opts = detectImportOptions(txtPath, 'FileType','text');
        % Forcer délimiteurs: tab + espace
        opts.Delimiter = {"\t"," "};
        opts = setvartype(opts, 'double');
        A = readmatrix(txtPath, opts);

        if isempty(A) || all(isnan(A), 'all')
            % Fallback
            error('EmptyOrNaN');
        end
        return
    catch
        % Fallback “manuel”
        txt  = fileread(txtPath);
        % Enlever guillemets au cas où:
        txt  = strrep(txt, '"', '');
        % Unifier tabs -> espace
        txt  = regexprep(txt, '[\t]+', ' ');
        % Réduire espaces multiples
        txt  = regexprep(txt, ' +', ' ');

        % Essai avec textscan
        fid = fopen(txtPath, 'r');
        C = textscan(fid, '%f', 'Delimiter', {' ','\t'}, 'MultipleDelimsAsOne', true);
        fclose(fid);
        v = C{1};
        if isempty(v)
            % Essai avec sscanf
            v = sscanf(txt, '%f');
        end
        % On suppose 4 ou 8 colonnes; si le fichier alterne, readmatrix gère mieux.
        % On reconstruit en colonnes en devinant le nombre de colonnes via la première ligne.
        firstLine = regexp(txt, '^[^\n\r]+', 'match', 'once', 'lineanchors');
        nCols = numel(strsplit(strtrim(firstLine)));
        if ~ismember(nCols, [4 8])
            % à défaut, essaie 8 puis 4
            if mod(numel(v), 8) == 0
                nCols = 8;
            elseif mod(numel(v), 4) == 0
                nCols = 4;
            else
                error('Impossible de déduire le nombre de colonnes pour %s', txtPath);
            end
        end
        A = reshape(v, nCols, []).';
    end
end
