#ifndef ARRAY_INCLUDED
#define ARRAY_INCLUDED

#define MISHA_DEBUG 0
//#define ASSERT( x ) { if( !( x ) ) _asm{ int 0x03 } }
  #define ASSERT( x ) { if( !( x ) ) __debugbreak(); }

#if MISHA_DEBUG
#include "Array.inl"
#define      Pointer( ... )      Array< __VA_ARGS__ >
#define ConstPointer( ... ) ConstArray< __VA_ARGS__ >
template< class C > void        VFreeArray( Array< C >& a )                                     { a.Free( ); }
template< class C > void         FreeArray( Array< C >& a )                                     { a.Free( ); }
template< class C > Array< C > VAllocArray( int size , int alignment=1 , const char* name=NULL ){ return Array< C >( size , alignment , false , true  , name ); }
template< class C > Array< C >  AllocArray( int size , int alignment=1 , const char* name=NULL ){ return Array< C >( size , alignment , false , false , name ); }
template< class C > Array< C > NullPointer( void )                                              { return Array< C >( ); }
template< class C > C*      PointerAddress( Array< C >& a )                                     { return a.pointer(); }
template< class C > Array< C >  GetPointer( C& c )                                              { return Array< C >::FromPointer( &c , 1 ); }
#else // !MISHA_DEBUG
#define      Pointer( ... )       __VA_ARGS__*
#define ConstPointer( ... ) const __VA_ARGS__*
#ifdef WIN32
template< class C > C*     AllocArray( int size , int alignment=1 , const char* name=NULL ){ return (C*)_aligned_malloc( sizeof(C) * size , alignment ); }
#define                     FreeArray( ptr ) if( ptr ) _aligned_free( ptr ) , ptr = NULL
#else // !WIN32
template< class C > C*     AllocArray( int size , int alignment=1 , const char* name=NULL ){ return (C*)memalign( sizeof(C) * size , alignment ); }
#define                     FreeArray( ptr ) if( ptr ) free( ptr ) , ptr = NULL
#endif // WIN32
#define                    VFreeArray( ptr ) if( ptr ) VirtualFree( ptr , 0 , MEM_RELEASE ) , ptr = NULL
template< class C > C*    VAllocArray( int size , int alignment=1 , const char* name=NULL ){ return (C*)VirtualAlloc( NULL , size*sizeof(C) , MEM_RESERVE | MEM_COMMIT , PAGE_READWRITE ); }
template< class C > C*    NullPointer( void ){ return NULL; }
template< class C > C* PointerAddress( C* c ){ return c; }
template< class C > C*     GetPointer( C& c ){ return &c; }
#endif // MISHA_DEBUG
#endif // ARRAY_INCLUDED
