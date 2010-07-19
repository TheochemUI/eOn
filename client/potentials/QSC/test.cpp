#include <stdio.h>
#include <string.h>
#include <iostream>

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

#include "QSC.h"

void init(void);
void report_errors(lua_State *L, int status);
int lua_set_r(lua_State *L);
int lua_set_z(lua_State *L);
int lua_set_box(lua_State *L);
int lua_get_force(lua_State *L);

QSC potential;
double *R;
long N;
long *atomicNrs;
double *box;
double *F;
double U;

int lua_get_force(lua_State *L)
{
    potential.force(N, R, atomicNrs, F, &U, box);
    lua_pushnumber(L, U);
    lua_newtable(L);
    for (int i=0; i<3*N; i++) {
        lua_pushnumber(L, i+1);
        lua_pushnumber(L, F[i]);
        lua_settable(L, -3);
    }
    return 2;
}

int lua_set_box(lua_State *L)
{
    luaL_checktype(L, 1, LUA_TTABLE);

    int n = luaL_getn(L, 1);
    box = new double[n];

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 1, i);
        box[i-1] = lua_tonumber(L, -1);
    }

    return 0;
}

int lua_set_z(lua_State *L)
{
    luaL_checktype(L, 1, LUA_TTABLE);

    int n = luaL_getn(L, 1);
    atomicNrs = new long[n];
    N = n;
    F = new double[3*n];

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 1, i);
        atomicNrs[i-1] = (long) lua_tonumber(L, -1);
    }

    return 0;
}

int lua_set_r(lua_State *L)
{
    luaL_checktype(L, 1, LUA_TTABLE);

    int n = luaL_getn(L, 1);
    R = new double[n];

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 1, i);
        double num = lua_tonumber(L, -1);
        R[i-1] = num;
    }

    return 0;
}


/*
void initialize(void) {
    //long N=4;
    //double R[] = { 50.0,50.0,50.0, 
    //               52.0,50.0,50.0,
    //               50.0,52.0,50.0,
    //               50.0,50.0,52.0};
    //const long atomicNrs[] = {46,79,46,79}; 
    //long N=2;
    //double R[] = { .0,.0,.0,
    //              3.0,.0,.0};
    //const long atomicNrs[] = {46,46};

    for (int i=0;i<75;i++) {
        pot.force(N, R, atomicNrs, F, &U, box);
        R[3] += .5;
        R[7] += .5;
        R[11] += .5;
        pot.force(N, R, atomicNrs, F, &U, box);
        printf("dR=%10.4g U=%10.4g\n", R[3]-R[0], U);
        printf("mag forces:\n");
        for (int j=0;j<N;j++) {
            double mag;
            mag = sqrt(F[3*j]*F[3*j] + F[3*j+1]*F[3*j+1] + 
                       F[3*j+2]*F[3*j+2]);
            printf("%13.4g\n", mag);
        }
        //for (int i=0;i<3*N;i++)
        //{
        //    if (i%3==0) printf("\n");
        //    printf("%10.4g ", F[i]); 
        //}
        //printf("\n");
    }
    pot.cleanMemory(); 
}
*/

void report_errors(lua_State *L, int status)
{
  if ( status!=0 ) {
    std::cerr << "-- " << lua_tostring(L, -1) << std::endl;
    lua_pop(L, 1); // remove error message
  }
}

int main(int argc, char **argv)
{
    for (int n=1; n<argc; ++n) {
        const char* file = argv[n];

        lua_State *L = lua_open();
        luaL_openlibs(L);

        lua_register(L, "set_r", lua_set_r);
        lua_register(L, "set_z", lua_set_z);
        lua_register(L, "set_box", lua_set_box);
        lua_register(L, "get_force", lua_get_force);

        //std::cerr << "-- Loading file: " << file << std::endl;

        int s = luaL_loadfile(L, file);

        if (s==0) {
            s = lua_pcall(L, 0, LUA_MULTRET, 0);
        }

        report_errors(L, s);

        lua_close(L);
    }


    return 0;
}
    
