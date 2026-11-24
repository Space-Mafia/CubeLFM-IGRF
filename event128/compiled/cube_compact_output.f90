module cube_compact_output
 use iso_fortran_env, only: real32, real64, int32
 use mpi_f08
 use lfm_cube_evolution, only: cube_model
 implicit none
 private
 public :: cube_output
 type cube_output
   character(:),allocatable :: prefix,stage,tool,python
   logical :: compact=.true.,first=.true.
 contains
   procedure :: open=>io_open, frame=>io_frame, checkpoint=>io_checkpoint, restore=>io_restore, close=>io_close
 end type
contains
 function quoted(value) result(word)
   character(*),intent(in)::value
   character(:),allocatable::word
   integer::i
   word="'"
   do i=1,len_trim(value)
     if(value(i:i)=="'")then
       word=word//"'"//achar(92)//"''"
     else
       word=word//value(i:i)
     endif
   enddo
   word=word//"'"
 end function
 function number(i) result(word)
   integer,intent(in)::i
   character(:),allocatable::word
   character(32)::s
   write(s,'(i0)')i
   word=trim(s)
 end function
 subroutine checked(code,what)
   integer,intent(in)::code
   character(*),intent(in)::what
   integer::st
   if(code/=0)then
     write(*,'(a,1x,i0)') 'CubeLFM IO: '//what,code
     call MPI_Abort(MPI_COMM_WORLD,1,st)
   endif
 end subroutine
 subroutine invoke(io,c,action,args)
   class(cube_output),intent(in)::io
   type(cube_model),intent(in)::c
   character(*),intent(in)::action
   character(*),intent(in),optional::args
   character(:),allocatable::cmd
   integer::code,status,st
   code=0
   if(c%rank==0)then
     cmd=quoted(io%python)//' '//quoted(io%tool)//' '//action//' '//quoted(io%prefix)
     if(present(args))cmd=cmd//' '//args
     call execute_command_line(cmd,wait=.true.,exitstat=code,cmdstat=status)
     if(status/=0)code=status
   endif
   call MPI_Bcast(code,1,MPI_INTEGER,0,MPI_COMM_WORLD,st)
   call checked(st,'broadcast IO status')
   call checked(code,action)
 end subroutine
 subroutine io_open(io,c,prefix)
   class(cube_output),intent(inout)::io
   type(cube_model),intent(in)::c
   character(*),intent(in)::prefix
   character(4096)::s
   integer::st,k
   if(c%rank==0)then
     call get_environment_variable('CUBELFM_OUTPUT_FORMAT',s,status=st)
     if(st==0.and.trim(s)/='legacy'.and.trim(s)/='compact')call checked(1,'unknown CUBELFM_OUTPUT_FORMAT')
     io%compact=(st/=0.or.trim(s)/='legacy')
   endif
   call MPI_Bcast(io%compact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,st)
   call checked(st,'format')
   io%prefix=trim(prefix)
   io%stage=io%prefix//'.io-tmp'
   if(.not.io%compact)return
   call get_environment_variable('CUBELFM_IO_PYTHON',s,status=st)
   if(st/=0.or.len_trim(s)==0)s='python3'
   io%python=trim(s)
   call get_environment_variable('CUBELFM_IO_TOOL',s,status=st)
   if(st/=0.or.len_trim(s)==0)then
     call get_command_argument(0,s,status=st)
     call checked(st,'executable path')
     k=index(trim(s),'/',back=.true.)
     s=s(:k)//'compact_io.py'
   endif
   io%tool=trim(s)
   call invoke(io,c,'init')
 end subroutine
 subroutine tile_stream(c,file,grid)
   type(cube_model),intent(in)::c
   character(*),intent(in)::file
   logical,intent(in)::grid
   type(MPI_File)::fh
   integer(MPI_OFFSET_KIND)::pos,bytes
   integer(int32)::header(7)
   integer::st,n(3),channel
   real(real32),allocatable::view(:,:,:)
   n=c%n
   header=int([n,c%layout%face,c%layout%first],int32)
   bytes=44_MPI_OFFSET_KIND
   if(grid)then
     bytes=bytes+8*size(c%p%x,kind=MPI_OFFSET_KIND)
   else
     bytes=bytes+32_MPI_OFFSET_KIND*product(n)+8*size(c%phi,kind=MPI_OFFSET_KIND)
   endif
   call MPI_Exscan(bytes,pos,1,MPI_OFFSET,MPI_SUM,MPI_COMM_WORLD,st)
   call checked(st,'rank byte offset')
   if(c%rank==0)pos=0
   call MPI_File_open(MPI_COMM_WORLD,file,ior(MPI_MODE_CREATE,ior(MPI_MODE_EXCL,MPI_MODE_WRONLY)),MPI_INFO_NULL,fh,st)
   call checked(st,'open shared staging file')
   call MPI_File_write_at_all(fh,pos,header,7,MPI_INTEGER4,MPI_STATUS_IGNORE,st)
   call checked(st,'tile layout')
   pos=pos+28
   call MPI_File_write_at_all(fh,pos,[c%time,c%rp],2,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,st)
   call checked(st,'tile clock')
   pos=pos+16
   if(grid)then
     call MPI_File_write_at_all(fh,pos,c%p%x,size(c%p%x),MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,st)
     call checked(st,'static coordinates')
   else
     allocate(view(n(1),n(2),n(3)))
     do channel=1,8
       if(channel<=5)then
         view=real(c%q(1:n(1),1:n(2),1:n(3),channel),real32)
       else
         view=real(c%q(1:n(1),1:n(2),1:n(3),channel)+c%background_cell(:,:,:,channel-5),real32)
       endif
       call MPI_File_write_at_all(fh,pos,view,size(view),MPI_REAL4,MPI_STATUS_IGNORE,st)
       call checked(st,'primitive field')
       pos=pos+4*size(view,kind=MPI_OFFSET_KIND)
     enddo
     call MPI_File_write_at_all(fh,pos,c%phi,size(c%phi),MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,st)
     call checked(st,'CT face flux')
   endif
   call MPI_File_close(fh,st)
   call checked(st,'close shared staging file')
 end subroutine
 subroutine io_frame(io,c)
   class(cube_output),intent(inout)::io
   type(cube_model),intent(inout)::c
   integer::st
   if(.not.io%compact)then
     call c%dump(io%prefix)
     call c%dump_mi(io%prefix)
     return
   endif
   if(io%first)call tile_stream(c,io%stage//'/grid',.true.)
   call tile_stream(c,io%stage//'/fields',.false.)
   call c%dump_mi(io%stage//'/raw')
   call MPI_Barrier(MPI_COMM_WORLD,st)
   call checked(st,'publish closed MHD and MI')
   call invoke(io,c,'frame',number(c%step)//' '//number(c%ranks)//' '//number(merge(1,0,c%mi_enabled)))
   io%first=.false.
 end subroutine
 subroutine io_checkpoint(io,c,fingerprint)
   class(cube_output),intent(in)::io
   type(cube_model),intent(inout)::c
   character(*),intent(in)::fingerprint
   if(.not.io%compact)then
     call c%checkpoint(io%prefix,fingerprint)
     return
   endif
   call c%checkpoint(io%stage//'/raw',fingerprint)
   call invoke(io,c,'checkpoint',number(c%step))
 end subroutine
 subroutine io_restore(io,c,path,fingerprint)
   class(cube_output),intent(in)::io
   type(cube_model),intent(inout)::c
   character(*),intent(in)::path,fingerprint
   integer::st,k
   k=len_trim(path)
   if(k>=3)then
     if(path(k-2:k)=='.h5')then
       if(.not.io%compact)call checked(1,'compact checkpoint requires compact output mode')
       call invoke(io,c,'restore',quoted(path))
       call c%restore(io%stage//'/restart',fingerprint)
       call MPI_Barrier(MPI_COMM_WORLD,st)
       call checked(st,'all ranks finished native restore')
       call invoke(io,c,'clean-restore')
       return
     endif
   endif
   call c%restore(path,fingerprint)
 end subroutine
 subroutine io_close(io,c)
   class(cube_output),intent(in)::io
   type(cube_model),intent(in)::c
   if(io%compact)call invoke(io,c,'end')
 end subroutine
end module
