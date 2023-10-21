from django.shortcuts import render, redirect
from django.http import HttpResponse
# from .forms import LoginForm, SignUpForm, ForgotPass, VerifyCodeForm, NewPassForm
from django.contrib import messages
# Create your views here.

def index(request):
    return render(request, 'pages/index.html',{"title": "CENTRO EDUCATIVO DE NIVEL SECUNDARIO N° 451"})