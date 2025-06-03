function mono_next_grlex_test ( )

%*****************************************************************************80
%
%% mono_next_grlex_test tests MONO_NEXT_GRLEX.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 February 2015
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'mono_next_grlex_test\n' );
  fprintf ( 1, '  MONO_NEXT_GRLEX returns the next monomial\n' );
  fprintf ( 1, '  in graded lexicographic order.\n' );

  m = 4;

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Let M =  %d\n', m );

  k = 0;
  x = zeros ( m, 1 );

  while ( true )

    d = sum ( x(1:m) );
    fprintf ( 1, '  %2d  %2d  |  %2d  %2d  %2d  %2d\n', k, d, x(1:m) );
    if ( x(1) == 3 )
      break
    end
    k = k + 1;
    x = mono_next_grlex ( m, x );
  end

  return
end
