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
int lua_get_new_sim(lua_State *L);
int lua_test(lua_State *L);

struct simulation {
    QSC *potential;
    double *R;
    long N;
    long *atomicNrs;
    double *box;
    double *F;
    double U;
};

int lua_get_new_sim(lua_State *L)
{
    luaL_checktype(L, 1, LUA_TNUMBER);
    long n = (long) lua_tonumber(L, 1);
    simulation *s = new simulation;
    s->N = n;
    s->potential = new QSC;
    s->R = new double[3*n];
    s->box = new double[3];
    s->atomicNrs = new long[n];
    s->F = new double[3*n];
    lua_pushlightuserdata(L, s);
    return 1;
}

int lua_get_force(lua_State *L)
{
    simulation *s;
    luaL_checktype(L, 1, LUA_TLIGHTUSERDATA);
    s = (simulation*) lua_topointer(L, 1);

    s->potential->force(s->N, s->R, s->atomicNrs, s->F, &(s->U), s->box);
    lua_pushnumber(L, s->U);
    lua_newtable(L);
    for (int i=0; i<3*(s->N); i++) {
        lua_pushnumber(L, s->F[i]);
        lua_rawseti(L, -2, i+1);
    }
    return 2;
}

int lua_set_box(lua_State *L)
{
    simulation *s;
    luaL_checktype(L, 1, LUA_TLIGHTUSERDATA);
    s = (simulation*) lua_topointer(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    int n = luaL_getn(L, 2);

    if (n != 3) {
        luaL_error(L, "set_box: box must be 3 dimensions");
    }

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 2, i);
        s->box[i-1] = lua_tonumber(L, -1);
    }

    return 0;
}

int lua_set_z(lua_State *L)
{
    simulation *s;
    luaL_checktype(L, 1, LUA_TLIGHTUSERDATA);
    s = (simulation*) lua_topointer(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    int n = luaL_getn(L, 2);

    if (n != s->N) {
        luaL_error(L, "set_z: must specify types for all atoms");
    }

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 2, i);
        s->atomicNrs[i-1] = (long) lua_tonumber(L, -1);
    }

    return 0;
}

int lua_set_r(lua_State *L)
{
    simulation *s;
    luaL_checktype(L, 1, LUA_TLIGHTUSERDATA);
    s = (simulation*) lua_topointer(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    int n = luaL_getn(L, 2);

    if (n != 3*s->N) {
        printf("%i %li\n", n, 3*s->N);
        return luaL_error(L, "set_r: table of incorrect length");
    }

    for (int i=1; i<=n; i++) {
        lua_rawgeti(L, 2, i);
        double num = lua_tonumber(L, -1);
        s->R[i-1] = num;
    }

    return 0;
}

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

        lua_register(L, "get_new_sim", lua_get_new_sim);
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
