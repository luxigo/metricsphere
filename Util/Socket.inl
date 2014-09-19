/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#define DEBUG_SOCKET 0


template<class C>
bool ReceiveOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Receiving (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( !ReceiveOnSocket( s ,  data , dataSize ) ) return false;
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Acknowledging (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack = 1;
		if( !SendOnSocket( s , &ack , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to acknowledge receive (%d)\n" , s );
			return false;
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
	return true;
}
template<class C>
bool SendOnSocket( SOCKET& s , ConstPointer( C ) data , int dataSize , bool blockingSend )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Sending (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , &ack , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement(%d)\n" , s );
			return false;
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
	if( !SendOnSocket( s , data , dataSize ) ) return false;
	return true;
}
template<class C>
bool SendOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Sending (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , &ack , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement(%d)\n" , s );
			return false;
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
	if( !SendOnSocket( s , data , dataSize ) ) return false;
	return true;
}


template<class C>
void ReceiveOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend , const char* errorMessage , ... )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Receiving (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( !ReceiveOnSocket( s ,  data , dataSize ) )
	{
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , errorMessage );
		vfprintf( stderr , errorMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
		exit(0);
	}
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Acknowledging (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack = 1;
		if( !SendOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to acknowledge receive (%d)\n" , s );
			{
				MyWinSock::StderrLock lock;
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
}

template<class C>
void SendOnSocket( SOCKET& s , ConstPointer( C ) data , int dataSize , bool blockingSend , const char* errorMessage , ... )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Sending (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement(%d)\n" , s );
			{
				MyWinSock::StderrLock lock;
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
	if( !SendOnSocket( s , data , dataSize ) )
	{
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , errorMessage );
		vfprintf( stderr , errorMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
		exit(0);
	}
}
template<class C>
void SendOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend , const char* errorMessage , ... )
{
#if DEBUG_SOCKET
	if( blockingSend ) printfId( " Sending (%d)\n" , s );
#endif // DEBUG_SOCKET
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement(%d)\n" , s );
			{
				MyWinSock::StderrLock lock;
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s );
#endif // DEBUG_SOCKET
	}
	if( !SendOnSocket( s , data , dataSize ) )
	{
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , errorMessage );
		vfprintf( stderr , errorMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
		exit(0);
	}
}

template<class C>
bool ReceiveOnSocket( SOCKET& s , Pointer( C ) data , int dataSize )
{
	int rec=0;
	while( rec!=dataSize )
	{
		int tmp = recv( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec , 0 );
		packetsReceived++;
		if (tmp <= 0)
		{
			MyWinSock::StderrLock lock;
			if(!tmp)	fprintfId( stderr , "Connection Closed...\n" );
			else		fprintfId( stderr , "recv from client failed (%d): %s\n" , s , LastSocketError() );
#if STORE_CONNECTION_TABLE
			if( s!=INVALID_SOCKET )
			{
				MyWinSock::SystemLock lock;
				if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
					fprintf( stderr , "Could not find connection value for %d\n" , s );
				else
				{
					ConnectionData& data = MyWinSock::ConnectionTable[ s ];
					fprintf( stderr , "\tConnection: %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
				}
			}
#endif // STORE_CONNECTION_TABLE
			return false;
		}
		rec+=tmp;
	}
	return true;
}
template<class C>
bool SendOnSocket( SOCKET& s , ConstPointer( C ) data , int dataSize )
{
	int sent = send( s , ( ConstPointer( char ) )data, dataSize, 0 );
	if (sent == SOCKET_ERROR)
	{
		MyWinSock::StderrLock lock;
		fprintfId( stderr , "send to client failed (%d): %s\n" , s , LastSocketError());
#if STORE_CONNECTION_TABLE
		if( s!=INVALID_SOCKET )
		{
			MyWinSock::SystemLock lock;
			if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
				fprintf( stderr , "Could not find connection value for %d\n" , s );
			else
			{
				ConnectionData& data = MyWinSock::ConnectionTable[ s ];
				fprintf( stderr , "\tConnection: %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
			}
		}
#endif // STORE_CONNECTION_TABLE
		return false;
	}
	packetsSent++;
	return true;
}
template<class C>
bool SendOnSocket( SOCKET& s , Pointer( C ) data , int dataSize )
{
	int sent = send( s , ( ConstPointer( char ) )data, dataSize, 0 );
	if (sent == SOCKET_ERROR)
	{
		MyWinSock::StderrLock lock;
		fprintfId( stderr , "send to client failed (%d): %s\n" , s , LastSocketError());
#if STORE_CONNECTION_TABLE
		if( s!=INVALID_SOCKET )
		{
			MyWinSock::SystemLock lock;
			if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
				fprintf( stderr , "Could not find connection value for %d\n" , s );
			else
			{
				ConnectionData& data = MyWinSock::ConnectionTable[ s ];
				fprintf( stderr , "\tConnection: %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
			}
		}
#endif // STORE_CONNECTION_TABLE
		return false;
	}
	packetsSent++;
	return true;
}
template<class C>
void ReceiveOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , const char* errorMessage , ... )
{
	int rec=0;
	while( rec!=dataSize )
	{
		int tmp = recv( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec , 0 );
		packetsReceived++;
		if (tmp <= 0)
		{
			MyWinSock::StderrLock lock;
			if(!tmp)	fprintfId( stderr , "Connection Closed...\n" );
			else		fprintfId( stderr , "recv from client failed (%d): %s\n" , s , LastSocketError() );
#if STORE_CONNECTION_TABLE
			if( s!=INVALID_SOCKET )
			{
				MyWinSock::SystemLock lock;
				if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
					fprintf( stderr , "Could not find connection value for %d\n" , s );
				else
				{
					ConnectionData& data = MyWinSock::ConnectionTable[ s ];
					fprintf( stderr , "\tConnection %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
				}
			}
#endif // STORE_CONNECTION_TABLE
			{
				fprintf( stderr , "\t" );
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
		rec+=tmp;
	}
}
template<class C>
void SendOnSocket( SOCKET& s , ConstPointer( C ) data , int dataSize , const char* errorMessage , ... )
{
	int sent = send( s , ( Pointer( char ) )data, dataSize, 0 );
	if (sent == SOCKET_ERROR)
	{
		MyWinSock::StderrLock lock;
		fprintfId( stderr , "send to client failed (%d): %s\n" , s , LastSocketError());
#if STORE_CONNECTION_TABLE
		if( s!=INVALID_SOCKET )
		{
			MyWinSock::SystemLock lock;
			if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
				fprintf( stderr , "Could not find connection value for %d\n" , s );
			else
			{
				ConnectionData& data = MyWinSock::ConnectionTable[ s ];
				fprintf( stderr , "\tConnection: %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
			}
		}
#endif // STORE_CONNECTION_TABLE
		{
			fprintf( stderr , "\t" );
			va_list args;
			va_start( args , errorMessage );
			vfprintf( stderr , errorMessage , args );
			va_end( args );
			fprintf( stderr , "\n" );
		}
		exit(0);
	}
	packetsSent++;
}
template< class C >
void SendOnSocket( SOCKET& s , Pointer( C ) data , int dataSize , const char* errorMessage , ... )
{
	int sent = send( s , ( Pointer( char ) )data, dataSize, 0 );
	if (sent == SOCKET_ERROR)
	{
		MyWinSock::StderrLock lock;
		fprintfId( stderr , "send to client failed (%d): %s\n" , s , LastSocketError());
#if STORE_CONNECTION_TABLE
		if( s!=INVALID_SOCKET )
		{
			MyWinSock::SystemLock lock;
			if( MyWinSock::ConnectionTable.find( s ) == MyWinSock::ConnectionTable.end( ) )
				fprintf( stderr , "Could not find connection value for %d\n" , s );
			else
			{
				ConnectionData& data = MyWinSock::ConnectionTable[ s ];
				fprintf( stderr , "\tConnection: %s:%d <-> %s:%d\n" , inet_ntoa( data.localAddr ) , data.localPort , inet_ntoa( data.peerAddr ) , data.peerPort );
			}
		}
#endif // STORE_CONNECTION_TABLE
		{
			fprintf( stderr , "\t" );
			va_list args;
			va_start( args , errorMessage );
			vfprintf( stderr , errorMessage , args );
			va_end( args );
			fprintf( stderr , "\n" );
		}
		exit(0);
	}
	packetsSent++;
}