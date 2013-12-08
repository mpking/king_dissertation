% Created: 08-19-2011
% Author: Michael Clemens - But don't blame me if it doesn't work! Many things on here are
% borrowed from other code.
%
% The VUPhys class is for producing phd dissertations 
% in the Vanderbilt department of physics and astronomy.  You can supply 
% the following optional arguments in the square brackets to 
% specify the thesis type:
%
%   phd     : Produces the PhD dissertation preliminary pages
%
% The default format is appropriate for printing, with blank pages 
% inserted after the preliminary pages in twoside mode so you can 
% send it directly to a two-sided printer. However, for ETD 
% submission the blank pages need to be removed from the final output.
% The following option does this:
%
%   etd     : Produces an electronic copy with no blank pages in the preliminary section
%
% The rest of the class options are the same as the regular book class.  
% A few to remember:
%
%   oneside : Produces single sided print layout (recommended for theses less than 50 pages)
%   twoside : Produces single sided print layout (the default if you remove oneside)
%
% The VUPhys class provides the following macros:
%
%   \makepreliminarypages : Makes the preliminary pages
%   \clearemptydoublepage : same as \cleardoublepage but doesn't put page numbers
%                           on blank intervening pages
%   \singlespace          : switch to single spaced lines
%   \doublespace          : switch to double spaced lines
%
% ------------------------------------------------------------------------------------------------------
%
\NeedsTeXFormat{LaTeX2e} \ProvidesClass{Dissertation}

% ---------------------------- declarations -------------------------
%
% These macros are used to declare arguments needed for the
% construction of the preliminary pages  

% The year and month the degree is awarded
\newcommand{\Year}[1]{\gdef\@Year{#1}}
\newcommand{\Month}[1]{\gdef\@Month{#1}}

% The full name of the degree
\newcommand{\degree}[1]{\gdef\@degree{#1}}

% The name of this document (thesis/dissertation)
\newcommand{\docname}[1]{\gdef\@docname{#1}}

% First line of title
\newcommand{\TitleTop}[1]{\gdef\@TitleTop{\mbox{\uppercase{#1}}}}

% Second line of title
\newcommand{\TitleBottom}[1]{\gdef\@TitleBottom{\mbox{\uppercase{#1}}}}

% Dedication text
\newcommand{\Dedication}[1]{\gdef\@Dedication{#1}}

% Acknowledgments text
\newcommand{\Acknowledgments}[1]{\gdef\@Acknowledgments{#1}}

% The author's name
\newcommand{\Author}[1]{\gdef\@Author{#1}}

% The name of the advisor
\newcommand{\Advisor}[1]{\gdef\@Advisor{#1}}

% The name of the committee member 2
\newcommand{\MemberA}[1]{\gdef\@MemberA{#1}}

% The name of the committee member 3
\newcommand{\MemberB}[1]{\gdef\@MemberB{#1}}

% The name of the committee member 4
\newcommand{\MemberC}[1]{\gdef\@MemberC{#1}}

% The name of the committee member 5
\newcommand{\MemberD}[1]{\gdef\@MemberD{#1}}

% Set default values for fields
  \Year{1905}
  \Month{January}
  \Author{Author}
  \TitleTop{First line of title}
  \TitleBottom{ } % default is empty
  \Acknowledgments{Acknowledgment text goes here.}
  \degree{DOCTOR OF PHILOSOPHY}
  \docname{Dissertation}
  \Advisor{Advisor}
  \MemberA{Committee Member A}
  \MemberB{Committee Member B}
  \MemberC{Committee Member C}
  \MemberD{Committee Member D}
  
%--------------------Superscript and Subscript Macros--------------

\newcommand{\super}[1]{\ensuremath{^{\textrm{#1}}}}
\newcommand{\sub}[1]{\ensuremath{_{\textrm{#1}}}}  
    
% ---------------------------- options ------------------------------

% A command to switch to single spaced lines
\newcommand{\singlespace}{\renewcommand{\baselinestretch}{1}\small\normalsize}

% A command to switch to double spaced lines
\newcommand{\doublespace}{\renewcommand{\baselinestretch}{1.66}\small\normalsize}

% A command pirated from chngpage.sty
\DeclareRobustCommand{\ch@ngetext}{%
  \setlength{\@colht}{\textheight}\setlength{\@colroom}{\textheight}%
  \setlength{\vsize}{\textheight}\setlength{\columnwidth}{\textwidth}%
  \if@twocolumn%
    \advance\columnwidth-\columnsep \divide\columnwidth\tw@%
    \@firstcolumntrue%
  \fi%
  \setlength{\hsize}{\columnwidth}%
  \setlength{\linewidth}{\hsize}%
}

% A command to make margins right for the initial single sided business.
\newcommand{\preliminarymargins}{%
    \addtolength{\textwidth}{-0.5in}%
    \addtolength{\evensidemargin}{0.5in}%
    \ch@ngetext%
    }

% A command to fix the margins after the initial single sided business.
\newcommand{\fixmargins}{%
    \addtolength{\textwidth}{0.5in}
    \addtolength{\evensidemargin}{-0.5in}
    \ch@ngetext%
}

% Changes to PhD preliminary section if the 'phd' option is specified
\DeclareOption{phd}{
  \degree{DOCTOR OF PHILOSOPHY}
  \docname{Dissertation}
  \newcommand{\makepreliminarypages}{ 
    \preliminarymargins
    \titlepage 
    \copyrightpage
    \dedicationpage 
    \acknowledgmentspage 
    \fixmargins
    \renewcommand{\clearemptydoublepage}{\cle@remptydoublep@ge}
  }
}

% --------------------- Some commands to handle the single sided preliminary pages ------------------

% Define the '\clearemptydoublepage' command to clear pages but not number any blank pages inserted.
% This is taken from the VUThesis.cls file
\let\cle@rdoublep@ge\cleardoublepage
\newcommand{\cle@remptydoublep@ge}{
  \clearpage
  \if@twoside
  \ifodd\c@page\else
  \fi\fi
  {\pagestyle{empty}\cle@rdoublep@ge}}
\newcommand{\clearemptydoublepage}{\cle@remptydoublep@ge}

% Redefine above commands if etd option is specified.  The blank pages make printing nice,
% but they don't want them in the submitted PDF
\DeclareOption{etd}{
    \renewcommand{\clearemptydoublepage}{ \clearpage }
  }

% ------------------------ Load the class and needed packages ---------------------------------

% Load the book class  
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}} 
\ProcessOptions\relax 
\LoadClass[letterpaper,12pt]{book}

% The afterpage package is required to make single sided formal pages 
% in a double sided environment
\RequirePackage{afterpage}

% Note: the hyperref package is required to make an appropriate ETD.
% However, we don't require it here since it is supposed to be the last
% package loaded and students may want to load other packages in the 
% main tex file.  So that this class file doesn't crash if the student
% forgets to load hyperref, we have used the following commands below:
%
%   \providecommand\phantomsection{}
%   \providecommand\pdfbookmark[3][]{}
%
% These commands provide dummy versions of the macros, but won't
% bother the real versions if the hyperref package is loaded in the 
% tex file.


% ---------------------------- main code ----------------------------

% Set the margins to VU specifications for the single sided preliminary
% pages.  At the end of the \makepreliminarypages command we fix the margins
% to be appropriate alternating values for double sided printing (if selected).
% If the \makepreliminarypages macro is not run, this never gets fixed. 
  \setlength{\marginparwidth}{0pt}
  \setlength{\marginparsep}{0pt}
  \setlength{\oddsidemargin}{0.5in}
  \setlength{\evensidemargin}{0in}
  \setlength{\textwidth}{6in}
  \setlength{\topmargin}{0in}
  \setlength{\headheight}{0.125in}
  \setlength{\headsep}{0.25in}
  \setlength{\textheight}{8.625in}
  \setlength{\footskip}{0.25in}
  \raggedbottom

% Redefine the Table of Contents to deal with some blank page 
% and bookmarking issues relating to ETD submission
\let\TEMPtableofcontents\tableofcontents
\renewcommand{\tableofcontents}{
  \clearemptydoublepage
  \providecommand\phantomsection{} \phantomsection
  %\addcontentsline{toc}{chapter}{Table of Contents}
  \TEMPtableofcontents
}
 
% Redefine the List of Figures to deal with some blank page 
% and bookmarking issues
\let\TEMPlistoffigures\listoffigures
\renewcommand{\listoffigures}{
  \providecommand\phantomsection{} \phantomsection
  \addcontentsline{toc}{chapter}{List of Figures}
  \TEMPlistoffigures
}

% Redefine the List of Tables to deal with some blank page 
% and bookmarking issues
\let\TEMPlistoftables\listoftables
\renewcommand{\listoftables}{
  \providecommand\phantomsection{} \phantomsection
  \addcontentsline{toc}{chapter}{List of Tables}
  \TEMPlistoftables
}

% Redefine the Bibliography to deal with a bookmarking issues
\let\TEMPbibliography\bibliography
\renewcommand{\bibliography}{
  \providecommand\phantomsection{} \phantomsection
  \addcontentsline{toc}{chapter}{Bibliography}
  \TEMPbibliography
}

%---------------------------- The Preliminary Page Definitions --------------------------
 
% ============================== Title Page ===============================
\renewcommand{\titlepage}{
    \thispagestyle{empty}
    \begin{center}
    \@TitleTop\\[\baselineskip]
    \@TitleBottom\\[\baselineskip]
    By\\[\baselineskip]
    \@Author\\[\baselineskip]
    Dissertation\\[\baselineskip]
    Submitted to the Faculty of the\\[\baselineskip]
    Graduate School of Vanderbilt University\\[\baselineskip]
    in partial fulfillment of the requirements\\[\baselineskip]
    for the degree of\\[\baselineskip]
    \@degree \\[\baselineskip]
    in\\[\baselineskip]
    Electrical Engineering\\[\baselineskip]
    \@Month~\@Year \\[\baselineskip]
    Nashville, Tennessee\\[\baselineskip]
    \vfill
    Approved:
    \vfill
    Professor \@Advisor
    \vfill
    Professor \@MemberA
    \vfill
    Professor \@MemberB
    \vfill   
    Professor \@MemberC
    \vfill
    Professor \@MemberD
    \end{center}
    \clearemptydoublepage
  }


% ======================== Copyright page ===============================
\newcommand{\copyrightpage}{
    \thispagestyle{empty}
    \providecommand\pdfbookmark[3][]{} \pdfbookmark[0]{Copyright}{bm:Copyright}    
    \vspace*{\fill}
    \begin{center}
    Copyright \copyright\ \@Year~\@Author\\[\baselineskip]
    All Rights Reserved\\[5\baselineskip]
    \end{center}
    \vspace{1in}
    \vfill
    \clearemptydoublepage
  }
  
 % ======================== Dedication page ===============================
\newcommand{\dedicationpage}{
    \thispagestyle{plain}
    \providecommand\pdfbookmark[3][]{} \pdfbookmark[0]{Dedication}{bm:Dedication}    
    \vspace*{\fill}
    \begin{center}
    \@Dedication
    \end{center}
    \vspace{1in}
    \vfill
    \clearemptydoublepage
  }

% ========================= Acknowledgments ==============================
\newcommand{\acknowledgmentspage}{
    \renewcommand{\baselinestretch}{1}\small\normalsize
    \thispagestyle{plain}
    \begin{center}
    \providecommand\pdfbookmark[3][]{} \pdfbookmark[0]{Acknowledgments}{bm:Acknowledge}    
    \vspace*{0.375in}
    ACKNOWLEDGMENTS\\[3\baselineskip]
    \end{center}
    
    \renewcommand{\baselinestretch}{1.66} \small\normalsize%
    \@Acknowledgments  
    \clearemptydoublepage
  }

