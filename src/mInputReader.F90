!!>
!!         programa de elementos finitos em fortran 90
!!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!!
!!         Pilato Jr , Vinicius Pessoa  {pilato@deepsoft.com.br, pessoa@deepsoft.com.br}
!!
!!         Desenvolvido por DeepSoft para LNCC/MCT
!!         Rio de Janeiro, 11.2013-04.2014
!!
!!=================================================================================

!> Modulo responsavel por reunir subrotinas para leitura do arquivo de entrada.
module mInputReader

    !> Armazena as linhas do arquivo de input.
    character(len=200), allocatable :: file_lines(:)
    !> Armazena o numero de linhas no arquivo.
    integer*4 number_of_lines

    contains

    !---------------------------------------------------------------------------------

    !> Le arquivo de input e armazena seu conteudo em um array.
    !! @param file_name Nome do arquivo a ser lido.
    subroutine readInputFileDS()
        use mIO,   only: iin

        implicit none

        integer*4 success, lines_count
        character(len=200) file_line

        integer*4 :: main_number_of_lines, main_number_of_includes, i, include_index, inc_nlines, inc_inc, merge_lines
        character(len=200), allocatable :: main_file_lines(:)
        character(len=200), allocatable :: include_files(:)
        character(len=200) include_file
        integer*4, allocatable :: include_indexes(:), include_number_of_lines(:)

        main_number_of_lines = 0
        main_number_of_includes = 0

        call analyzeFileInput(main_number_of_lines, main_number_of_includes)

        if (main_number_of_includes.eq.0) then
            call createSimpleInputFile()
            return
        end if

        allocate(main_file_lines(main_number_of_lines))
        allocate(include_indexes(main_number_of_includes))
        allocate(include_number_of_lines(main_number_of_includes))
        allocate(include_files(main_number_of_includes))

        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            main_file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)

        number_of_lines = main_number_of_lines

        !Number of lines
        do i=1, main_number_of_includes
            include_index = findInclude(i, main_file_lines, main_number_of_lines)
            read(main_file_lines(include_index), '(A)') include_file

            include_file = adJustl(include_file)
            call analyzeFile(include_file, inc_nlines, inc_inc)
            number_of_lines = number_of_lines + inc_nlines
            include_indexes(i) = include_index
            include_number_of_lines(i) = inc_nlines
            include_files(i) = include_file
        end do

        !Prepare final struct.
        call prepareFileLines(include_indexes, include_number_of_lines, main_number_of_includes, main_file_lines)

        !Merge contensts.
        merge_lines = 0
        do i=1, main_number_of_includes
            call mergeIncludeContents(include_files(i), include_indexes(i) + merge_lines)
            merge_lines = merge_lines + include_number_of_lines(i)
        end do

        deallocate(main_file_lines)
        deallocate(include_indexes)
        deallocate(include_number_of_lines)
        return
    end subroutine readInputFileDS !***********************************************************************************


    !> Cria a estrutura de input usando um arquivo de entrada sem includes
    !! @param file_name Nome do arquivo a ser lido.
    subroutine createSimpleInputFile()
        use mIO,   only: iin

        implicit none

        integer*4 success, lines_count
        character(len=200) file_line

        number_of_lines = 0

        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
        end do
        rewind(iin)

        allocate(file_lines(number_of_lines))
        !TO-DO avoid two-times read
        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)
    end subroutine createSimpleInputFile !*****************************************************************************

    !> Le o conteudo do arquivo de include e armazena no array principal.
    !! @param   include_index   O index do include.
    !! @param   include_files   Array com includes.
    !! @param   include_line    A linha do include.
    subroutine mergeIncludeContents(include_file, include_line)

        implicit none

        integer*4 include_line
        character(len=200) include_file

        character(len=200) file_line
        integer*4 file_channel, success, current_index

        file_channel = 1

        current_index = include_line

        open(unit=file_channel, file=include_file)
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(current_index) = file_line
            current_index = current_index + 1
        end do
        close(file_channel)
    end subroutine mergeIncludeContents !******************************************************************************

    !> Efetua a alocacao da estrutura definitiva, preparando a linha dos arquivos originais para receber os includes
    !! @param   include_indexes             Array os indices de ocorrencias dos includes.
    !! @param   include_number_of_lines     Array com o numero de linhas de cada include
    !! @param   number_of_includes          Numero de includes.
    !! @param   original_file_lines         Linhas do arquivo de entrada original.
    subroutine prepareFileLines(include_indexes, include_number_of_lines, number_of_includes, original_file_lines)
        integer*4 number_of_includes, number_of_original_lines, line_index, shift_lines
        integer*4 include_indexes(:), include_number_of_lines(:)
        character(len=200) original_file_lines(:)

        integer*4 current_include_index, original_index

        allocate(file_lines(number_of_lines))

        current_include_index = 1
        original_index = 1
        line_index = 1
        shift_lines = 0
        do while ( line_index <= number_of_lines)
            if (original_index.eq.(include_indexes(current_include_index))) then
                line_index = line_index + include_number_of_lines(current_include_index)
                current_include_index = current_include_index + 1
            end if
            file_lines(line_index) = original_file_lines(original_index)
            line_index = line_index + 1
            original_index = original_index + 1
        end do


    end subroutine prepareFileLines !**********************************************************************************

    !> Efetua algumas analises no arquivo recebido.
    !! @param   number_of_lines     Numero de linhas.
    !! @param   number_of_include   Numero de ocorrencias da palavra include.
    subroutine analyzeFileInput(number_of_lines, number_of_includes)
        use mIO,   only: iin

        character(len=200) file_line
        integer*4 number_of_lines, number_of_includes

        character(len=50) include_keyword, formated_keyword
        integer*4 keyword_len, success

        include_keyword = "include"
        keyword_len = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines = 0
        number_of_includes = 0

!        lunitInicial = 15
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        rewind(iin)

    end subroutine analyzeFileInput !***************************************************************************************

    !> Efetua algumas analises no arquivo recebido.
    !! @param   file_name           O nome do arquivo.
    !! @param   number_of_lines     Numero de linhas.
    !! @param   number_of_include   Numero de ocorrencias da palavra include.
    subroutine analyzeFile(file_name, number_of_lines, number_of_includes)
        character(len=200) file_name, file_line
        integer*4 number_of_lines, number_of_includes

        character(len=50) include_keyword, formated_keyword
        integer*4 keyword_len, file_channel, success

        include_keyword = "include"
        keyword_len = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines = 0
        number_of_includes = 0

        file_channel = 2
        lunitInicial = 15
        file_channel = lunitInicial

        open(unit=file_channel, file=file_name)
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        close(file_channel)

    end subroutine analyzeFile !***************************************************************************************

    !> Procura a n-esima palavra-chave include.
    !! @param  position         Corresponde a posicao desejada.
    !! @param  file_lines       Linhas do arquivo.
    !! @param  number_of_lines  Numero de linhas atuais.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findInclude(position, file_lines, number_of_lines)
        implicit none
        integer*4 position, number_of_lines, current_position
        character(len=200) file_lines(:)
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len

        keyword = "include"
        keyword_len = len(trim(keyword)) + 2
        formated_keyword = trim('*' // trim(keyword) // '{')
        current_position = 0

        do i=1, number_of_lines
            file_line = file_lines(i)
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                current_position = current_position + 1
                if (current_position.eq.position) then
                    findInclude = i + 1
                    return
                end if
            end if
        end do
        findInclude = 0
        return
    end function findInclude !*****************************************************************************************

    !> Procura uma palavra-chave.
    !! @param  keyword A palavra-chave.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findKeyword(keyword)
        implicit none
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len
        do i=1, number_of_lines, 1
            file_line = file_lines(i)
            keyword_len = len(trim(keyword)) + 2
            formated_keyword = trim('*' // trim(keyword) // '{')
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                findKeyword = i + 1
                return
            end if
        end do
        findKeyword = 0
        return
    end function findKeyword !*****************************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo inteiro. Se nao encontrado, associa o valor default fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o valor inteiro sera atribuido.
    !! @param default_value Valor default.
    subroutine readIntegerKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        integer*4 target, default_value, keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) target
        return
    end subroutine readIntegerKeywordValue !***************************************************************************

    !> Efetua a leitura de uma palavra-chave do tipo array de inteiro. Se nao encontrado, associa o valor default fornecido.
    !! Obs.: Atentar para o fato dessa sub-rotina ter um do "infinito".
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o valor inteiro sera atribuido.
    !! @param default_value Valor default.
    !! @author Diego Volpatto
    subroutine readIntArrayKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        integer*4 default_value, keyword_line, lines_num, i
        integer*4, allocatable :: target(:)
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            lines_num = 1
            allocate(target(lines_num));
            target = default_value
            return
        end if
        lines_num = 0
        do
            file_line = adjustL(trim(file_lines(keyword_line+lines_num)))
            if (file_line .eq. '}') then
                if (lines_num .ge. 1) allocate(target(lines_num));
                exit
            else
                lines_num = lines_num + 1
            endif
        enddo
        do i=1,lines_num
            file_line = adjustL(trim(file_lines(keyword_line+i-1)))
            read(file_line, *) target(i)
        enddo
        return
    end subroutine readIntArrayKeywordValue !***************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo string. Se nao encontrado, associa o valor default fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde a string sera atribuido.
    !! @param default_value Valor default.
    subroutine readStringKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        character(len=*)          :: target, default_value
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        read(file_lines(keyword_line), '(a)') target
        return
    end subroutine readStringKeywordValue !****************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo real. Se nao encontrado, associa o valor default fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o real sera atribuido.
    !! @param default_value Valor default.
    subroutine readRealKeywordValue(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        real(8) target, default_value
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            target = default_value
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) target
        return
    end subroutine readRealKeywordValue !****************************************************************************

    !> Efetua a leitura de uma palavra-chave do tipo de um array bidimensional real. 
    !! A leitura eh realizada linha por linha.
    !! Se nao encontrado, associa o valor default fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde os valores serao atribuido.
    !! @param default_value Valor default.
    !! @author Diego T. Volpatto
    subroutine readRealMatrixValues(keyword, target, default_value)
        implicit none
        character(50) keyword
        character(120) file_line
        integer :: idx, i, j, n, m
        real(8)  default_value
        real*8, allocatable :: target(:,:)
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            allocate(target(1,1))
            target = default_value
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) n, m
        allocate(target(n,m))
        do i=1,n
            file_line = adjustL(trim(file_lines(keyword_line+i)))
            read(file_line, *) idx, (target(idx,j), j=1,m)
        enddo
        return
    end subroutine  !****************************************************************************

    !> Efetua a leitura de uma palavra-chave do tipo de um array bidimensional real. 
    !! A leitura eh realizada linha por linha.
    !! Se nao encontrado, associa o valor default fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param kbc           Tipo da CC (1 = Dirichlet, 2 = Neumann).
    !! @param vbc           Valor prescrito na CC.
    !! @param default_value Valor default.
    !! @author Diego T. Volpatto
    subroutine readBoundaryConditions(keyword, kbc, vbc, default_value)
        use mIO, only: iout
        implicit none
        character(50) keyword
        character(120) file_line
        integer :: idx, i, j, n, m
        real(8)  default_value
        integer*4, allocatable :: kbc(:)
        real*8, allocatable :: vbc(:)
        integer*4 keyword_line
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            allocate(kbc(1)); allocate(vbc(1))
            kbc = 0; vbc = 0.0d0
            return
        end if
        file_line = adjustL(trim(file_lines(keyword_line)))
        read(file_line, *) n
        allocate(kbc(n)); kbc = -1
        allocate(vbc(n)); vbc = 0.0d0
        write(iout,*) "Boundary values data:"
        do i=1,n
            file_line = adjustL(trim(file_lines(keyword_line+i)))
            read(file_line, *) idx, kbc(idx), vbc(idx)
            write(iout, *) "Boundary:", idx, kbc(idx), vbc(idx)
        enddo
        write(iout,*)
        return
    end subroutine  !****************************************************************************

end module mInputReader

